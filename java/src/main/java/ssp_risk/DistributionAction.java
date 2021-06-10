package ssp_risk;

@FunctionalInterface
interface DistributionAction {
  void accept(int key, double probability);
}
