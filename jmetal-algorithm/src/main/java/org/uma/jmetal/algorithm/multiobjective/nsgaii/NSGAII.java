package org.uma.jmetal.algorithm.multiobjective.nsgaii;

import org.uma.jmetal.algorithm.impl.AbstractGeneticAlgorithm;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.operator.impl.selection.RankingAndCrowdingSelection;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.SolutionListUtils;
import org.uma.jmetal.util.comparator.DominanceComparator;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.solution.DoubleSolution;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 */
@SuppressWarnings("serial")
public class NSGAII<S extends Solution<?>> extends AbstractGeneticAlgorithm<S, List<S>> {
  protected final int maxEvaluations;

  protected final SolutionListEvaluator<S> evaluator;

  protected int evaluations;
  protected Comparator<S> dominanceComparator ;

  protected int matingPoolSize;
  protected int offspringPopulationSize ;

  /**
   * Constructor
   */
  public NSGAII(Problem<S> problem, int maxEvaluations, int populationSize,
                int matingPoolSize, int offspringPopulationSize,
      CrossoverOperator<S> crossoverOperator, MutationOperator<S> mutationOperator,
      SelectionOperator<List<S>, S> selectionOperator, SolutionListEvaluator<S> evaluator) {
    this(problem, maxEvaluations, populationSize, matingPoolSize, offspringPopulationSize,
            crossoverOperator, mutationOperator, selectionOperator, new DominanceComparator<S>(), evaluator);
  }
  /**
   * Constructor
   */
  public NSGAII(Problem<S> problem, int maxEvaluations, int populationSize,
                int matingPoolSize, int offspringPopulationSize,
                CrossoverOperator<S> crossoverOperator, MutationOperator<S> mutationOperator,
      SelectionOperator<List<S>, S> selectionOperator, Comparator<S> dominanceComparator,
                SolutionListEvaluator<S> evaluator) {
    super(problem);
    this.maxEvaluations = maxEvaluations;
    setMaxPopulationSize(populationSize); ;

    this.crossoverOperator = crossoverOperator;
    this.mutationOperator = mutationOperator;
    this.selectionOperator = selectionOperator;

    this.evaluator = evaluator;
    this.dominanceComparator = dominanceComparator ;

    this.matingPoolSize = matingPoolSize ;
    this.offspringPopulationSize = offspringPopulationSize ;
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
    evaluations = getMaxPopulationSize();
  }

  @Override protected void updateProgress() {
    evaluations += offspringPopulationSize ;
  }

  @Override protected boolean isStoppingConditionReached() {
    return evaluations >= maxEvaluations;
  }

  @Override protected List<S> evaluatePopulation(List<S> population) {
    population = evaluator.evaluate(population, getProblem());

    return population;
  }

  /**
   * This method iteratively applies a {@link SelectionOperator} to the population to fill the mating pool population.
   *
   * @param population
   * @return The mating pool population
   */
  @Override
  protected List<S> selection(List<S> population) {
    List<S> matingPopulation = new ArrayList<>(population.size());
    for (int i = 0; i < matingPoolSize; i++) {
      S solution = selectionOperator.execute(population);
      matingPopulation.add(solution);
    }

    return matingPopulation;
  }

  /**
   * This methods iteratively applies a {@link CrossoverOperator} a  {@link MutationOperator} to the population to
   * create the offspring population. The population size must be divisible by the number of parents required
   * by the {@link CrossoverOperator}; this way, the needed parents are taken sequentially from the population.
   *
   * The number of solutions returned by the {@link CrossoverOperator} must be equal to the offspringPopulationSize
   * state variable
   *
   * @param matingPool
   * @return The new created offspring population
   */
  @Override
  protected List<S> reproduction(List<S> matingPool) {
    int numberOfParents = crossoverOperator.getNumberOfRequiredParents() ;

    checkNumberOfParents(matingPool, numberOfParents);

    List<S> offspringPopulation = new ArrayList<>(offspringPopulationSize);
    for (int i = 0; i < matingPool.size(); i += numberOfParents) {
      List<S> parents = new ArrayList<>(numberOfParents);
      for (int j = 0; j < numberOfParents; j++) {
        parents.add(population.get(i+j));
      }

      List<S> offspring = crossoverOperator.execute(parents);

      for(S s: offspring){
        mutationOperator.execute(s);
        offspringPopulation.add(s);
        if (offspringPopulation.size() >= offspringPopulationSize)
          break;
      }
    }
    return offspringPopulation;
  }

  @Override protected List<S> replacement(List<S> population, List<S> offspringPopulation) {
    List<S> jointPopulation = new ArrayList<>();
    jointPopulation.addAll(population);
    jointPopulation.addAll(offspringPopulation);

    RankingAndCrowdingSelection<S> rankingAndCrowdingSelection ;
    rankingAndCrowdingSelection = new RankingAndCrowdingSelection<S>(getMaxPopulationSize(), dominanceComparator) ;

    return rankingAndCrowdingSelection.execute(jointPopulation) ;
  }

  @Override public List<S> getResult() {
    return SolutionListUtils.getNondominatedSolutions(getPopulation());
  }

  @Override public String getName() {
    return "NSGAII" ;
  }

  @Override public String getDescription() {
    return "Nondominated Sorting Genetic Algorithm version II" ;
  }
}
