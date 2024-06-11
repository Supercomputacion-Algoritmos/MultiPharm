package org.uma.jmetal.runner.multiobjective;

import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.wasfga.WASFGA;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.operator.impl.crossover.PMXCrossover;
import org.uma.jmetal.operator.impl.mutation.PermutationSwapMutation;
import org.uma.jmetal.operator.impl.selection.BinaryTournamentSelection;
import org.uma.jmetal.problem.PermutationProblem;
import org.uma.jmetal.problem.multiobjective.MultiobjectiveTSP;
import org.uma.jmetal.solution.PermutationSolution;
import org.uma.jmetal.util.AbstractAlgorithmRunner;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.comparator.RankingAndCrowdingDistanceComparator;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;
import org.uma.jmetal.problem.multiobjective.MultiPharm;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.operator.impl.crossover.SBXCrossover;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.operator.impl.selection.BinaryTournamentSelection;
import org.uma.jmetal.util.ProblemUtils;
import org.uma.jmetal.operator.impl.crossover.DifferentialEvolutionCrossover;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class WASFGARunner extends AbstractAlgorithmRunner {
  /**
   * @param args Command line arguments.
   * @throws JMetalException
   * @throws FileNotFoundException
   * Invoking command:
  java org.uma.jmetal.runner.multiobjective.WASFGABinaryRunner problemName [referenceFront]
   */
  public static void main(String[] args) throws JMetalException, IOException {
    Problem<DoubleSolution> problem;
    //PermutationProblem<PermutationSolution<Integer>> problem;
    Algorithm<List<DoubleSolution>> algorithm;
    CrossoverOperator<DoubleSolution> crossover;
    MutationOperator<DoubleSolution> mutation;
    SelectionOperator<List<DoubleSolution>, DoubleSolution> selection;
    /*Algorithm<List<PermutationSolution<Integer>>> algorithm;
    PermutationProblem<PermutationSolution<Integer>> problem;
    CrossoverOperator<PermutationSolution<Integer>> crossover;
    MutationOperator<PermutationSolution<Integer>> mutation;
    SelectionOperator<List<PermutationSolution<Integer>>, PermutationSolution<Integer>> selection;*/
    
    String problemName ;
    if (args.length == 3){
      problemName = args[0];
      MultiPharm.queryPath = args[1];
      MultiPharm.targetPath = args[2];
    } else {
      problemName = "org.uma.jmetal.problem.multiobjective.MultiPharm";
      //referenceParetoFront = "jmetal-problem/src/test/resources/pareto_fronts/ZDT1.pf";
    }

    //problem = new MultiobjectiveTSP("/tspInstances/kroA100.tsp", "/tspInstances/kroB100.tsp");
    problem = (DoubleProblem) ProblemUtils.<DoubleSolution>loadProblem(problemName);

    //crossover = new PMXCrossover(0.9) ;

    //double mutationProbability = 0.2 ;
    //mutation = new PermutationSwapMutation<Integer>(mutationProbability) ;

    //selection = new BinaryTournamentSelection<PermutationSolution<Integer>>(new RankingAndCrowdingDistanceComparator<PermutationSolution<Integer>>());
    String referenceParetoFront = "" ;
    List<Double> referencePoint = null;

    ////problem = ProblemUtils.<DoubleSolution> loadProblem(problemName);
    //problem = new MultiobjectiveTSP("/tspInstances/kroA100.tsp", "/tspInstances/kroB100.tsp");

    referencePoint = new ArrayList<>();
    referencePoint.add(-0.5);
    referencePoint.add(-0.5);

    double crossoverProbability = 0.9 ;
    double crossoverDistributionIndex = 20.0 ;
    crossover = new SBXCrossover(crossoverProbability, crossoverDistributionIndex) ;
    
    //double cr = 1.0 ;
    //double f = 0.5 ;
    //crossover = new DifferentialEvolutionCrossover(cr, f, "rand/1/bin");

    double mutationProbability = 1.0 / problem.getNumberOfVariables() ;
    double mutationDistributionIndex = 20.0 ;
    mutation = new PolynomialMutation(mutationProbability, mutationDistributionIndex) ;

    selection = new BinaryTournamentSelection<DoubleSolution>(new RankingAndCrowdingDistanceComparator<DoubleSolution>());

    double epsilon = 0.01 ;
    algorithm = new WASFGA<DoubleSolution>(problem, 300, 2000, crossover, mutation, selection,
            new SequentialSolutionListEvaluator<DoubleSolution>(),epsilon, referencePoint) ;

    AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm)
            .execute() ;

    List<DoubleSolution> population = algorithm.getResult() ;
    long computingTime = algorithmRunner.getComputingTime() ;

    JMetalLogger.logger.info("Total execution time: " + computingTime + "ms");

    printFinalSolutionSet(population, algorithm.getName());
    if (!referenceParetoFront.equals("")) {
      printQualityIndicators(population, referenceParetoFront) ;
    }
  }
}
