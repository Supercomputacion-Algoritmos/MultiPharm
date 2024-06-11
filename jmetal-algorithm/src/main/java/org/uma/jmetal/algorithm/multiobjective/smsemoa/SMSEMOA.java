package org.uma.jmetal.algorithm.multiobjective.smsemoa;

import org.uma.jmetal.algorithm.impl.AbstractGeneticAlgorithm;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.qualityindicator.impl.Hypervolume;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.solutionattribute.Ranking;
import org.uma.jmetal.util.solutionattribute.impl.DominanceRanking;
import org.uma.jmetal.solution.DoubleSolution;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 */
@SuppressWarnings("serial")
public class SMSEMOA<S extends Solution<?>> extends AbstractGeneticAlgorithm<S, List<S>> {
  protected final int maxEvaluations;
  protected final double offset ;

  protected int evaluations;

  private Hypervolume<S> hypervolume;
  protected Comparator<S> dominanceComparator ;

  /**
   * Constructor
   */
  public SMSEMOA(Problem<S> problem, int maxEvaluations, int populationSize, double offset,
      CrossoverOperator<S> crossoverOperator, MutationOperator<S> mutationOperator,
      SelectionOperator<List<S>, S> selectionOperator, Comparator<S> dominanceComparator, Hypervolume<S> hypervolumeImplementation) {
    super(problem) ;
    this.maxEvaluations = maxEvaluations;
    setMaxPopulationSize(populationSize);

    this.offset = offset ;

    this.crossoverOperator = crossoverOperator;
    this.mutationOperator = mutationOperator;
    this.selectionOperator = selectionOperator;
    this.dominanceComparator = dominanceComparator ;
    this.hypervolume = hypervolumeImplementation ;
  }
  
  // New method by Savins
      @Override protected List<S> createInitialPopulation() {
        List<S> population = new ArrayList<>(maxPopulationSize);
        double molecules_per_Axis = 80.0;
        double angle = 6.283185 / molecules_per_Axis;
        double angle_temp = 0.0;

        // starting position
        DoubleSolution newSolution = (DoubleSolution) problem.createSolution();
        newSolution.setVariableValue(0, 0.0);
        newSolution.setVariableValue(1, 1.0);
        newSolution.setVariableValue(2, 0.0);
        newSolution.setVariableValue(3, 0.0);
        newSolution.setVariableValue(4, 0.0);
        newSolution.setVariableValue(5, 0.0);
        newSolution.setVariableValue(6, 0.0);
        newSolution.setVariableValue(7, 0.0);
        newSolution.setVariableValue(8, 0.0);
        newSolution.setVariableValue(9, 0.0);
        problem.evaluate((S) newSolution);
        population.add((S) newSolution);

        // Poses rotate around X axis
        angle_temp = angle;
        for (int i = 0; i < molecules_per_Axis; i++) {
            newSolution = (DoubleSolution) problem.createSolution();
            newSolution.setVariableValue(0, angle_temp);
            newSolution.setVariableValue(1, 1.0);
            newSolution.setVariableValue(2, 0.0);
            newSolution.setVariableValue(3, 0.0);
            newSolution.setVariableValue(4, 0.0);
            newSolution.setVariableValue(5, 0.0);
            newSolution.setVariableValue(6, 0.0);
            newSolution.setVariableValue(7, 0.0);
            newSolution.setVariableValue(8, 0.0);
            newSolution.setVariableValue(9, 0.0);
            problem.evaluate((S) newSolution);
            population.add((S) newSolution);
            angle_temp += angle;
        }

        // Poses rotate around Y axis
        angle_temp = angle;
        for (int i = 0; i < molecules_per_Axis; i++) {
            newSolution = (DoubleSolution) problem.createSolution();
            newSolution.setVariableValue(0, angle_temp);
            newSolution.setVariableValue(1, 0.0);
            newSolution.setVariableValue(2, 1.0);
            newSolution.setVariableValue(3, 0.0);
            newSolution.setVariableValue(4, 0.0);
            newSolution.setVariableValue(5, 0.0);
            newSolution.setVariableValue(6, 0.0);
            newSolution.setVariableValue(7, 0.0);
            newSolution.setVariableValue(8, 0.0);
            newSolution.setVariableValue(9, 0.0);
            problem.evaluate((S) newSolution);
            population.add((S) newSolution);
            angle_temp += angle;
        }
        // Poses rotate around Z axis
        angle_temp = angle;
        for (int i = 0; i < molecules_per_Axis; i++) {
            newSolution = (DoubleSolution) problem.createSolution();
            newSolution.setVariableValue(0, angle_temp);
            newSolution.setVariableValue(1, 0.0);
            newSolution.setVariableValue(2, 0.0);
            newSolution.setVariableValue(3, 1.0);
            newSolution.setVariableValue(4, 0.0);
            newSolution.setVariableValue(5, 0.0);
            newSolution.setVariableValue(6, 0.0);
            newSolution.setVariableValue(7, 0.0);
            newSolution.setVariableValue(8, 0.0);
            newSolution.setVariableValue(9, 0.0);
            problem.evaluate((S) newSolution);
            population.add((S) newSolution);
            angle_temp += angle;
        }
        // Poses generated randomly

        for (int i = 0; i < maxPopulationSize-(1+molecules_per_Axis*3); i++) {
            newSolution = (DoubleSolution)problem.createSolution();
            problem.evaluate((S) newSolution);
            population.add((S) newSolution);
        }
        return population;
    }
    // END OF NEW METHOD BY SAVINS

  @Override protected void initProgress() {
    evaluations = getMaxPopulationSize() ;
  }

  @Override protected void updateProgress() {
    evaluations++ ;
  }

  @Override protected boolean isStoppingConditionReached() {
    return evaluations >= maxEvaluations ;
  }

  @Override protected List<S> evaluatePopulation(List<S> population) {
    for (S solution : population) {
      getProblem().evaluate(solution);
    }
    return population ;
  }

  @Override protected List<S> selection(List<S> population) {
    List<S> matingPopulation = new ArrayList<>(2);
    for (int i = 0; i < 2; i++) {
      S solution = selectionOperator.execute(population);
      matingPopulation.add(solution);
    }

    return matingPopulation;
  }

  @Override protected List<S> reproduction(List<S> population) {
    List<S> offspringPopulation = new ArrayList<>(1);

    List<S> parents = new ArrayList<>(2);
    parents.add(population.get(0));
    parents.add(population.get(1));

    List<S> offspring = crossoverOperator.execute(parents);

    mutationOperator.execute(offspring.get(0));

    offspringPopulation.add(offspring.get(0));
    return offspringPopulation;
  }

  @Override protected List<S> replacement(List<S> population, List<S> offspringPopulation) {
    List<S> jointPopulation = new ArrayList<>();
    jointPopulation.addAll(population);
    jointPopulation.addAll(offspringPopulation);

    Ranking<S> ranking = computeRanking(jointPopulation);
    List<S> lastSubfront = ranking.getSubfront(ranking.getNumberOfSubfronts()-1) ;

    lastSubfront = hypervolume.computeHypervolumeContribution(lastSubfront, jointPopulation) ;

    List<S> resultPopulation = new ArrayList<>() ;
    for (int i = 0; i < ranking.getNumberOfSubfronts()-1; i++) {
      for (S solution : ranking.getSubfront(i)) {
        resultPopulation.add(solution);
      }
    }

    for (int i = 0; i < lastSubfront.size()-1; i++) {
      resultPopulation.add(lastSubfront.get(i)) ;
    }

    return resultPopulation ;
  }

  @Override public List<S> getResult() {
    return getPopulation();
  }

  protected Ranking<S> computeRanking(List<S> solutionList) {
    Ranking<S> ranking = new DominanceRanking<S>(dominanceComparator);
    ranking.computeRanking(solutionList);

    return ranking;
  }

  @Override public String getName() {
    return "SMSEMOA" ;
  }

  @Override public String getDescription() {
    return "S metric selection EMOA" ;
  }
}