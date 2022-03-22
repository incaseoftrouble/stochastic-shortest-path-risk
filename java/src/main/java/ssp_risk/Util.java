package ssp_risk;

import java.util.function.IntToDoubleFunction;

public final class Util {
  private static final double EPS = 1.0e-07;
  private static final double LOOSE = 1.0e-6;

  private Util() {
  }

  public static boolean doublesEqual(double d1, double d2) {
    return doubleZero(d1 - d2);
  }

  public static boolean doublesEqualLoose(double d1, double d2) {
    return doubleZeroLoose(d1 - d2);
  }

  public static boolean doubleZero(double d) {
    return Math.abs(d) < EPS;
  }

  public static boolean doubleZeroLoose(double d) {
    return Math.abs(d) < LOOSE;
  }

  public static boolean doubleLessOrEqualZeroLoose(double d) {
    return d < LOOSE;
  }

  public static boolean doubleLessOrEqualZero(double d) {
    return d < EPS;
  }

  public static boolean doublesLessOrEqual(double d1, double d2) {
    return d1 - d2 <= EPS;
  }

  public static boolean doublesLessOrEqualLoose(double d1, double d2) {
    return d1 - d2 <= LOOSE;
  }

  public static double sum(int from, int to, IntToDoubleFunction values) {
    double sum = 0.0;
    double c = 0.0;
    for (int i = from; i < to; i++) {
      double y = values.applyAsDouble(i) - c;
      double t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
    return sum;
  }
}
