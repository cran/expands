package core.analysis.ngs.algorithms;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import core.utils.Common;

public class ExPANdS {

  public static final int DEFAULT_MAX_PM = 6;
  private static int min_PM = 1;;
  private static int PN = 2;

  /**
   * frequency of mutated cell
   */
  private List<FitAlternatives> solutions;

  /**
   * Observed/measured copy number
   */
  private double CN;

  /**
   * Observed/measured allele frequency
   */
  private double AF;

  /**
   * Estimated allele frequency -error
   */
  private double CN_Error = 0.1;

  /**
   * Observed/measured ploidy of B allele in non-utated cells
   */
  private int PN_B;

  private double minErr;

  private double[] alternativeCN;
  private FitAlternatives bestSolution;
  private int max_PM;

  /**
   * 
   * @return the deviation of the observed from the expected copy number and
   *         allele frequency, for the inferred solution
   */
  public double getDeviation() {
    return minErr;
  }

  public ExPANdS(double af, double cn, int pnb, int max_PM) {
    this.max_PM = max_PM;
    bestSolution = new FitAlternatives(1, -1, -1);
    bestSolution.add(0, Double.NaN, Double.NaN);
    this.AF = af;
    this.CN = cn;
    this.PN_B = pnb;
    double step = 0.005;
    int alt = (int) Math.floor(((CN + CN_Error) - Math.max(0, CN - CN_Error))
        / step);
    this.alternativeCN = new double[alt];
    alternativeCN[0] = Math.max(0, CN - CN_Error);
    for (int i = 1; i < alternativeCN.length; i++) {
      alternativeCN[i] = alternativeCN[i - 1] + step;
    }
  }
 

  public void run() {
    minErr = Double.POSITIVE_INFINITY;
    int nIter = max_PM - min_PM + 1;
    int[][] PM = new int[nIter][nIter];
    int[][] PM_B = new int[nIter][nIter];
    solutions = new ArrayList<ExPANdS.FitAlternatives>(PM.length * PM_B.length);
    for (int pm = min_PM; pm <= max_PM; pm++) {

      double copyPenality = Math.pow(Math.abs(pm - Math.round(CN)) + 1, 3);

      for (int pmb = min_PM; pmb <= pm; pmb++) {
        if (pmb / (double) pm <= PN_B / (double) PN) {
          continue;
        }
        RealMatrix coefficients = new Array2DRowRealMatrix(new double[][] {
            { pm - PN }, { pmb - PN_B } }, false);

        // Assign
        int i = pm - min_PM;
        int j = pmb - min_PM;
        PM[i][j] = pm;
        PM_B[i][j] = pmb;

        // FitAlternatives for same PM and PM_B considering error in CN
        // measurement
        FitAlternatives alternatives = new FitAlternatives(
            alternativeCN.length, pm, pmb);
        for (int k = 0; k < alternativeCN.length; k++) {
          double cn = alternativeCN[k];
          RealVector constants = new ArrayRealVector(new double[] { cn - PN,
              cn * AF - PN_B }, false);
          DecompositionSolver solver = new SingularValueDecomposition(
              coefficients).getSolver();
          RealVector solution = solver.solve(constants);
          RealVector errorterm = coefficients.operate(solution).subtract(
              constants);
          double sumErr = copyPenality
              * (Math.abs(errorterm.getEntry(0))
                  + Math.abs(errorterm.getEntry(1)) + Math.abs(CN - cn));
          alternatives.add(k, solution.getEntry(0), sumErr);
          if (solution.getEntry(0) > -0.1 && solution.getEntry(0) <= 1.1
              && minErr > sumErr) {
            minErr = sumErr;
            this.bestSolution = alternatives;
            //System.out.println(solution.getEntry(0)+"\t"+errorterm.getEntry(0)
            // +
            // "\t"+errorterm.getEntry(1));
          }
        }

        solutions.add(alternatives);

      }
    }
  }

  public int getPM_B() {
    return bestSolution.getPM_B();
  }

  public int getPM() {
    return bestSolution.getPM();
  }

  public double getF() {
    try {
      return bestSolution.getBest().getF();
    } catch (Exception e) {
      return Double.NaN;
    }
  }

  public void printFitsToFile(File f) throws IOException {
    BufferedWriter w = Common.getWriter(f.getAbsolutePath());
    w.write(Common.toString(SOLUTION_ENTITIES, "\t"));
    w.newLine();
    for (FitAlternatives fits : this.solutions) {
      Iterator<Solution> iter = fits.iterator();
      while (iter.hasNext()) {

        w.write(iter.next().toString());
        w.newLine();

      }
    }
    w.flush();
    w.close();
  }

  public Collection<FitAlternatives> solutions() {
    return solutions;
  }

  public static void main(String[] args) {

    File[] files = (new File(System.getProperty("user.dir"))
        .listFiles(new FilenameFilter() {

          @Override
          public boolean accept(File arg0, String arg1) {
            return arg1.endsWith(".snv");
          }
        }));
    Arrays.sort(files);
    for (File f : files) {
      System.out.println(f.getName());

      File pair = new File(f.getAbsolutePath() + ".expands");
      if (pair.exists()) {
        System.out.println("ExPANdS file for " + f.getName()
            + " already exists. Skipped.");
        continue;
      }

      int count = 0;
      try {
        BufferedReader r = Common.getReader(f.getAbsolutePath());
        BufferedWriter w = Common.getWriter(pair.getAbsolutePath());
        String[] header = r.readLine().trim().split("\\s+");
        w.write(Common.toString(header, "\t"));
        w.newLine();
        int countI = Common.firstIndexOf("Count", header);
        int afI = Common.firstIndexOf("AF_Tumor", header);
        int cnI = Common.firstIndexOf("CN_Estimate", header);
        int pnbI = Common.firstIndexOf("PN_B", header);
        int pABBBI = Common.firstIndexOf("pABBB", header);
        int pAABBI = Common.firstIndexOf("pAABB", header);
        int pAAABI = Common.firstIndexOf("pAAAB", header);
        int fI = Common.firstIndexOf("f", header);
        int pmI = Common.firstIndexOf("PM", header);
        int pmbI = Common.firstIndexOf("PM_B", header);
        int devI = Common.firstIndexOf("dev", header);
        for (String l = r.readLine(); l != null; l = r.readLine()) {
          String[] features = l.trim().split("\t");
          if (features.length > header.length) {
            features = Arrays.copyOfRange(features, 1, features.length);
          }
          int id = (int) Double.parseDouble(features[countI]);
          int pnb = (int) Math.round(Double.parseDouble(features[pABBBI]));
          try {
            ExPANdS expands = new ExPANdS(Double.parseDouble(features[afI]),
                Double.parseDouble(features[cnI]), pnb, DEFAULT_MAX_PM);
            expands.run();
            double pAABB = Double.parseDouble(features[pAABBI]);
            double pABBB = Double.parseDouble(features[pABBBI]);
            double pAAAB = Double.parseDouble(features[pAAABI]);
            if (pABBB >= 0.5 || pAABB + pAAAB >= 0.9) {
              expands.printFitsToFile(new File(f.getAbsolutePath() + "." + id
                  + ".fit"));
              count++;
            }
            features[pmI] = "" + expands.getPM();
            features[pmbI] = "" + expands.getPM_B();
            features[fI] = "" + expands.getF();
            features[devI] = "" + expands.getDeviation();
            features[pnbI] = "" + pnb;
          } catch (Exception e) {
            e.printStackTrace();
          }

          w.write(Common.toString(features, "\t").replace("Infinity", "Inf"));
          w.newLine();

        }
        w.flush();
        w.close();
      } catch (FileNotFoundException e1) {
        // TODO Auto-generated catch block
        e1.printStackTrace();
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }

    }
  }

  public class FitAlternatives implements Iterator<Solution>,
      Iterable<Solution> {

    private double[] f;
    private double[] dev;
    private int pm;
    private int pmb;
    int idx = -2;

    private FitAlternatives(int numerAlternatives, int pm, int pmb) {
      f = new double[numerAlternatives];
      dev = new double[numerAlternatives];
      this.pm = pm;
      this.pmb = pmb;
    }

    public int getPM_B() {
      return pmb;
    }

    public int getPM() {
      return pm;
    }

    public Iterator<Solution> iterator() {
      idx = -1;
      return this;
    }

    public double getF(int k) {
      return f[k];
    }

    public double getDev(int k) {
      return dev[k];
    }

    public int size() {
      return dev.length;
    }

    public Solution getBest() {
      int bestIdx = Common.argmin(dev);
      return new Solution(PN, PN_B, pm, pmb, CN, AF, f[bestIdx], dev[bestIdx]);
    }

    public void add(int idx, double f, double dev) {
      this.f[idx] = f;
      this.dev[idx] = dev;
    }

    @Override
    public boolean hasNext() {
      return idx > -2 && idx < f.length - 1;
    }

    @Override
    public Solution next() {
      idx++;
      if (idx >= dev.length) {
        return null;
      }
      return new Solution(PN, PN_B, pm, pmb, CN, AF, f[idx], dev[idx]);
    }

    @Override
    public void remove() {
      // TODO Auto-generated method stub

    }

  }

  public static final String[] SOLUTION_ENTITIES = { "PN", "PN_B", "PM",
      "PM_B", "CN_Estimate", "AF_Tumor", "f", "dev" };

  public final class Solution {

    private double f;
    private double AF;
    private double CN;
    private int PM_B;
    private int PM;
    private int PN_B;
    private int PN;
    private double dev;

    public Solution(int PN, int PN_B, int PM, int PM_B, double CN, double AF,
        double f, double dev) {
      this.PN = PN;
      this.PN_B = PN_B;
      this.PM = PM;
      this.PM_B = PM_B;
      this.CN = CN;
      this.AF = AF;
      this.f = f;
      this.dev = dev;

    }

    public double getF() {
      return f;
    }

    public double[] toDouble() {
      return new double[] { PN, PN_B, PM, PM_B, CN, AF, f, dev };
    }

    public String toString() {
      return Common.toString(toDouble(), "\t");
    }
  }
}
