package ssp_risk;

public final class Util {
  private static final double EPS = 1.0e-16;
  private static final double LOOSE = 1.0e-6;

  private Util() {
  }

  public static boolean doublesEqual(double d1, double d2) {
    return doubleZero(d1 - d2);
  }

  public static boolean doubleZero(double d) {
    return Math.abs(d) < EPS;
  }

  public static boolean lessZeroLoose(double d) {
    return d < LOOSE;
  }
}
