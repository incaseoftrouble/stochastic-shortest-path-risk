package ssp_risk;

public final class Util {
  private static final double EPS = 1.0e-14;
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
}
