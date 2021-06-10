package ssp_risk;

import com.google.gson.annotations.SerializedName;
import java.time.Duration;
import java.util.Map;
import javax.annotation.Nullable;

public record Result(
    InputData input,
    ModelData model,
    ProblemParameters problem,
    @SerializedName("ssp")
    SspSolution sspSolution,
    @Nullable
    @SerializedName("vi")
    ValueIterationSolution viSolution,
    @Nullable
    @SerializedName("lp")
    LinearProgrammingSolution lpSolution
) {
  public record InputData(
      String name,
      String[] args
  ) {}

  public record ModelData(
      int states,
      long actions,
      long transitions,
      @SerializedName("initial_state")
      int initialState,
      @SerializedName("absorbing_states")
      int absorbingStates
  ) {}

  public record ProblemParameters(
      @SerializedName("goal_states")
      int goalStates,
      double[] thresholds
  ) {}

  public record SspSolution(
      @SerializedName("construction_time")
      Duration constructionTime,
      @SerializedName("total_time")
      Duration totalTime,
      double result
  ) {}

  public record ThresholdResult(
      int var,
      double cvar,
      Duration time
  ) {}

  public record ValueIterationSolution(
      Duration time,
      Map<Double, ThresholdResult> results
  ) {}

  public record LinearProgrammingSolution(
      Duration time,
      @SerializedName("preprocessing_time")
      Duration reachPreprocessingTime,
      @SerializedName("construction_time")
      Duration constructionTime,
      @SerializedName("solving_time")
      Duration solvingTime,
      Map<Double, ThresholdResult> results
  ) {}
}
