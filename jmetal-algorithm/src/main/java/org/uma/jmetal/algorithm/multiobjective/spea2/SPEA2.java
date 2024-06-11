package org.uma.jmetal.algorithm.multiobjective.spea2;

import org.uma.jmetal.algorithm.impl.AbstractGeneticAlgorithm;
import org.uma.jmetal.algorithm.multiobjective.spea2.util.EnvironmentalSelection;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.solutionattribute.impl.StrengthRawFitness;
import org.uma.jmetal.solution.DoubleSolution;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Juan J. Durillo
 **/
@SuppressWarnings("serial")
public class SPEA2<S extends Solution<?>> extends AbstractGeneticAlgorithm<S, List<S>> {
  protected final int maxIterations;
  protected final SolutionListEvaluator<S> evaluator;
  protected int iterations;
  protected List<S> archive;
  protected final StrengthRawFitness<S> strenghtRawFitness = new StrengthRawFitness<S>();
  protected final EnvironmentalSelection<S> environmentalSelection;
  protected final int k ;

  public SPEA2(Problem<S> problem, int maxIterations, int populationSize,
      CrossoverOperator<S> crossoverOperator, MutationOperator<S> mutationOperator,
      SelectionOperator<List<S>, S> selectionOperator, SolutionListEvaluator<S> evaluator,
               int k) {
    super(problem);
    this.maxIterations = maxIterations;
    this.setMaxPopulationSize(populationSize);

    this.k = k ;
    this.crossoverOperator = crossoverOperator;
    this.mutationOperator = mutationOperator;
    this.selectionOperator = selectionOperator;
    this.environmentalSelection = new EnvironmentalSelection<S>(populationSize, k);

    this.archive = new ArrayList<>(populationSize);

    this.evaluator = evaluator;
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

  @Override
  protected void initProgress() {
    iterations = 1;
  }

  @Override
  protected void updateProgress() {
    iterations++;
  }

  @Override
  protected boolean isStoppingConditionReached() {
    return iterations >= maxIterations;
  }

  @Override
  protected List<S> evaluatePopulation(List<S> population) {
    population = evaluator.evaluate(population, getProblem());
    return population;
  }

  @Override
  protected List<S> selection(List<S> population) {
    List<S> union = new ArrayList<>(2*getMaxPopulationSize());
    union.addAll(archive);
    union.addAll(population);
    strenghtRawFitness.computeDensityEstimator(union);
    archive = environmentalSelection.execute(union);
    return archive;
  }

  @Override
  protected List<S> reproduction(List<S> population) {
    List<S> offSpringPopulation= new ArrayList<>(getMaxPopulationSize());

    while (offSpringPopulation.size() < getMaxPopulationSize()){
      List<S> parents = new ArrayList<>(2);
      S candidateFirstParent = selectionOperator.execute(population);
      parents.add(candidateFirstParent);
      S candidateSecondParent;
      candidateSecondParent = selectionOperator.execute(population);
      parents.add(candidateSecondParent);

      List<S> offspring = crossoverOperator.execute(parents);
      mutationOperator.execute(offspring.get(0));
      offSpringPopulation.add(offspring.get(0));
    }
    return offSpringPopulation;
  }

  @Override
  protected List<S> replacement(List<S> population,
      List<S> offspringPopulation) {
    return offspringPopulation;
  }

  @Override
  public List<S> getResult() {
    return archive;
  }

  @Override public String getName() {
    return "SPEA2" ;
  }

  @Override public String getDescription() {
    return "Strength Pareto. Evolutionary Algorithm" ;
  }
}
