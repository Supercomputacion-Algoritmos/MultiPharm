package org.uma.jmetal.algorithm.multiobjective.moead;

import org.uma.jmetal.algorithm.multiobjective.moead.util.MOEADUtils;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.impl.crossover.DifferentialEvolutionCrossover;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.DoubleSolution;

import java.util.ArrayList;
import java.util.List;

/**
 * Class implementing the MOEA/D-DE algorithm described in :
 * Hui Li; Qingfu Zhang, "Multiobjective Optimization Problems With Complicated Pareto Sets,
 * MOEA/D and NSGA-II," Evolutionary Computation, IEEE Transactions on , vol.13, no.2, pp.284,302,
 * April 2009. doi: 10.1109/TEVC.2008.925798
 *
 * @author Antonio J. Nebro
 * @version 1.0
 */
@SuppressWarnings("serial")
public class MOEADMiriam extends AbstractMOEAD<DoubleSolution> {
    protected DifferentialEvolutionCrossover differentialEvolutionCrossover;

    public MOEADMiriam(Problem<DoubleSolution> problem,
                 int populationSize,
                 int resultPopulationSize,
                 int maxEvaluations,
                 MutationOperator<DoubleSolution> mutation,
                 CrossoverOperator<DoubleSolution> crossover,
                 FunctionType functionType,
                 String dataDirectory,
                 double neighborhoodSelectionProbability,
                 int maximumNumberOfReplacedSolutions,
                 int neighborSize) {
        super(problem, populationSize, resultPopulationSize, maxEvaluations, crossover, mutation, functionType,
                dataDirectory, neighborhoodSelectionProbability, maximumNumberOfReplacedSolutions,
                neighborSize);

        differentialEvolutionCrossover = (DifferentialEvolutionCrossover) crossoverOperator;
    }

    @Override
    public void run() {
        initializePopulation();
        initializeUniformWeight();
        initializeNeighborhood();
        idealPoint.update(population);
        ;

        evaluations = populationSize;
        do {
            int[] permutation = new int[populationSize];
            MOEADUtils.randomPermutation(permutation, populationSize);

            for (int i = 0; i < populationSize; i++) {
                int subProblemId = permutation[i];

                NeighborType neighborType = chooseNeighborType();
                List<DoubleSolution> parents = parentSelection(subProblemId, neighborType);

                differentialEvolutionCrossover.setCurrentSolution(population.get(subProblemId));
                List<DoubleSolution> children = differentialEvolutionCrossover.execute(parents);

                DoubleSolution child = children.get(0);
                mutationOperator.execute(child);
                problem.evaluate(child);

                evaluations++;

                idealPoint.update(child.getObjectives());
                updateNeighborhood(child, subProblemId, neighborType);
            }
        } while (evaluations < maxEvaluations);

    }



/* Original initializePopulation Savins
  protected void initializePopulation() {
    population = new ArrayList<>(populationSize);
    for (int i = 0; i < populationSize; i++) {
      DoubleSolution newSolution = (DoubleSolution)problem.createSolution();

      problem.evaluate(newSolution);
      population.add(newSolution);
    }
  }
  */
    // MIRIAM
    protected void initializePopulation() {
        population = new ArrayList<>(populationSize);
        double molecules_per_Axis = 80.0;
        double angle = 6.283185 / molecules_per_Axis;
        double angle_temp = 0.0;

       
        DoubleSolution newSolution = (DoubleSolution) problem.createSolution();
         newSolution.setVariableValue(0, 0.0);
         newSolution.setVariableValue(1, 0.0);
         newSolution.setVariableValue(2, Math.PI/2);
         newSolution.setVariableValue(3, 0.0);
         newSolution.setVariableValue(4, 0.0);
         newSolution.setVariableValue(5, 0.0);
         problem.evaluate(newSolution);
         population.add(newSolution);

      // Poses rotate around X axis
         //angle_temp = angle;
         for (int i = 0; i < molecules_per_Axis; i++) {
        	 newSolution = (DoubleSolution) problem.createSolution();
        	 newSolution.setVariableValue(0, angle_temp);
             newSolution.setVariableValue(1, 0.0);
             newSolution.setVariableValue(2, Math.PI/2);
             newSolution.setVariableValue(3, 0.0);
             newSolution.setVariableValue(4, 0.0);
             newSolution.setVariableValue(5, 0.0);
             problem.evaluate(newSolution);
             population.add(newSolution);
             angle_temp += angle;
       
         }
         
      // Poses rotate around Y axis
         angle_temp = angle;
         
         for (int i = 0; i < molecules_per_Axis; i++) {
        	 newSolution = (DoubleSolution) problem.createSolution();
        	 newSolution.setVariableValue(0, angle_temp);
             newSolution.setVariableValue(1, Math.PI/2);
             newSolution.setVariableValue(2, Math.PI/2);
             newSolution.setVariableValue(3, 0.0);
             newSolution.setVariableValue(4, 0.0);
             newSolution.setVariableValue(5, 0.0);
             problem.evaluate(newSolution);
             population.add(newSolution);
             angle_temp += angle;
             
         }
         
         // Poses rotate around Z axis
         angle_temp = angle;
         for (int i = 0; i < molecules_per_Axis; i++) {
        	 newSolution = (DoubleSolution) problem.createSolution();
        	 newSolution.setVariableValue(0, angle_temp);
             newSolution.setVariableValue(1, 0.0);
             newSolution.setVariableValue(2, 0.0);
             newSolution.setVariableValue(3, 0.0);
             newSolution.setVariableValue(4, 0.0);
             newSolution.setVariableValue(5, 0.0);
             problem.evaluate(newSolution);
             population.add(newSolution);
             angle_temp += angle;
         	
             
         }
         
         // RANDOM
         for (int i = 0; i < populationSize-(1+molecules_per_Axis*3); i++) {
             newSolution = (DoubleSolution)problem.createSolution();
             problem.evaluate(newSolution);
             population.add(newSolution);
         }
         
        // Poses rotate around X axis
        /*angle_temp = angle;
        for (int i = 0; i < molecules_per_Axis; i++) {
            newSolution = (DoubleSolution) problem.createSolution();
            newSolution.setVariableValue(0, angle_temp);
            newSolution.setVariableValue(1, 1.0);
            newSolution.setVariableValue(2, 0.0);
            newSolution.setVariableValue(3, 0.0);
            newSolution.setVariableValue(4, 0.0);
            newSolution.setVariableValue(5, 0.0);
            //newSolution.setVariableValue(6, 0.0);
            //newSolution.setVariableValue(7, 0.0);
            //newSolution.setVariableValue(8, 0.0);
            //newSolution.setVariableValue(9, 0.0);
            problem.evaluate(newSolution);
            population.add(newSolution);
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
            //newSolution.setVariableValue(6, 0.0);
            //newSolution.setVariableValue(7, 0.0);
            //newSolution.setVariableValue(8, 0.0);
            //newSolution.setVariableValue(9, 0.0);
            problem.evaluate(newSolution);
            population.add(newSolution);
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
            //newSolution.setVariableValue(6, 0.0);
            //newSolution.setVariableValue(7, 0.0);
            //newSolution.setVariableValue(8, 0.0);
            //newSolution.setVariableValue(9, 0.0);
            problem.evaluate(newSolution);
            population.add(newSolution);
            angle_temp += angle;
        }
        // Poses generated randomly

        for (int i = 0; i < populationSize-(1+molecules_per_Axis*3); i++) {
            newSolution = (DoubleSolution)problem.createSolution();
            problem.evaluate(newSolution);
            population.add(newSolution);
        }*/
       /* DoubleSolution newSolution = (DoubleSolution) problem.createSolution();
        for (int i = 0; i < populationSize; i++) {
            newSolution = (DoubleSolution)problem.createSolution();
            problem.evaluate(newSolution);
            population.add(newSolution);
        }*/
    }


    @Override
    public String getName() {
        return "MOEADMiriam";
    }

    @Override
    public String getDescription() {
        return "Multi-Objective Evolutionary Algorithm based on Decomposition";
    }
}
