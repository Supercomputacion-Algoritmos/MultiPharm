package org.uma.jmetal.algorithm.multiobjective.wasfga;

import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.operator.impl.selection.BinaryTournamentSelection;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.AlgorithmBuilder;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;

import java.util.List;

/**
 * @author Miriam R. Ferrandez
 */
public class WASFGABuilder<S extends Solution<?>> implements AlgorithmBuilder<WASFGA<S>> {
  /**
   * WASFGABuilder class
   */
  protected final Problem<S> problem;
  protected int maxIterations;
  protected int populationSize;
  protected CrossoverOperator<S> crossoverOperator;
  protected MutationOperator<S> mutationOperator;
  protected SelectionOperator<List<S>, S> selectionOperator;
  protected SolutionListEvaluator<S> evaluator;
  protected double epsilon ;
  protected List<Double> referencePoint;

  /**
   * WASFGABuilder constructor
   */
  public WASFGABuilder(Problem<S> problem, CrossoverOperator<S> crossoverOperator,
      MutationOperator<S> mutationOperator, SelectionOperator<List<S>, S> selectionOperator, List<Double> referencePoint) {
    this.problem = problem;
    maxIterations = 500;
    populationSize = 300;
    this.crossoverOperator = crossoverOperator ;
    this.mutationOperator = mutationOperator ;
    this.selectionOperator = selectionOperator ;
    evaluator = new SequentialSolutionListEvaluator<S>();
    epsilon = 0.01 ;
    this.referencePoint = referencePoint;
  }

  public WASFGABuilder<S> setMaxIterations(int maxIterations) {
    if (maxIterations < 0) {
      throw new JMetalException("maxIterations is negative: " + maxIterations);
    }
    this.maxIterations = maxIterations;

    return this;
  }

  public WASFGABuilder<S> setPopulationSize(int populationSize) {
    if (populationSize < 0) {
      throw new JMetalException("Population size is negative: " + populationSize);
    }

    this.populationSize = populationSize;

    return this;
  }

  public WASFGABuilder<S> setSelectionOperator(SelectionOperator<List<S>, S> selectionOperator) {
    if (selectionOperator == null) {
      throw new JMetalException("selectionOperator is null");
    }
    this.selectionOperator = selectionOperator;

    return this;
  }

  public WASFGABuilder<S> setSolutionListEvaluator(SolutionListEvaluator<S> evaluator) {
    if (evaluator == null) {
      throw new JMetalException("evaluator is null");
    }
    this.evaluator = evaluator;

    return this;
  }

  public WASFGABuilder<S> setEpsilon(double epsilon) {
    this.epsilon = epsilon ;

    return this;
  }

  public WASFGA<S> build() {
    WASFGA<S> algorithm = null ;
    algorithm = new WASFGA<S>(problem, populationSize, maxIterations, crossoverOperator,
          mutationOperator, selectionOperator, evaluator, epsilon, referencePoint);
    
    return algorithm ;
  }

  /* Getters */
  public Problem<S> getProblem() {
    return problem;
  }

  public int getMaxIterations() {
    return maxIterations;
  }

  public int getPopulationSize() {
    return populationSize;
  }

  public CrossoverOperator<S> getCrossoverOperator() {
    return crossoverOperator;
  }

  public MutationOperator<S> getMutationOperator() {
    return mutationOperator;
  }

  public SelectionOperator<List<S>, S> getSelectionOperator() {
    return selectionOperator;
  }

  public SolutionListEvaluator<S> getSolutionListEvaluator() {
    return evaluator;
  }
}
