package org.uma.jmetal.runner.multiobjective;

import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.moead.AbstractMOEAD;
import org.uma.jmetal.algorithm.multiobjective.moead.MOEADBuilder;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.impl.crossover.DifferentialEvolutionCrossover;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.problem.multiobjective.MultiPharmMRF;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.AbstractAlgorithmRunner;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.ProblemUtils;

import java.io.FileNotFoundException;
import java.util.List;

/**
 * Class for configuring and running the MOEA/D algorithm
 *
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 */
public class MOEADRunnerMultiPharmMRF extends AbstractAlgorithmRunner {

  public static String queryName;
  public static String targetName;

  /**
   * @param args Command line arguments.
   * @throws SecurityException
   * Invoking command:
  java org.uma.jmetal.runner.multiobjective.MOEADRunner problemName [referenceFront]
   */
  public static void main(String[] args) throws FileNotFoundException {
    DoubleProblem problem;
    Algorithm<List<DoubleSolution>> algorithm;
    MutationOperator<DoubleSolution> mutation;
    DifferentialEvolutionCrossover crossover;

    String problemName ;
    String referenceParetoFront = "";
    if (args.length == 1) {
      problemName = args[0];
    } else if (args.length == 2) {
      problemName = args[0];
      referenceParetoFront = args[1];
    }else if (args.length == 3){

      problemName = args[0];
      MultiPharmMRF.queryPath = args[1];
      MultiPharmMRF.targetPath = args[2];
    } else if (args.length == 5){
      problemName = args[0];
      MultiPharmMRF.queryPath = args[1];
      MultiPharmMRF.targetPath = args[2];
      queryName  = args[3];
      targetName = args[4];
    } else {
      problemName = "org.uma.jmetal.problem.multiobjective.MultiPharmMRF";

      //referenceParetoFront = "/Users/savins/temp/pruebaJMETAL/jMetal-jmetal-5.8/jmetal-problem/src/test/resources/pareto_fronts/DTLZ2.3D.pf";
    }

    problem = (DoubleProblem)ProblemUtils.<DoubleSolution> loadProblem(problemName);

    double cr = 1.0 ;
    double f = 0.5 ;
    crossover = new DifferentialEvolutionCrossover(cr, f, "rand/1/bin");

    double mutationProbability = 1.0 / problem.getNumberOfVariables();
    double mutationDistributionIndex = 20.0;
    mutation = new PolynomialMutation(mutationProbability, mutationDistributionIndex);

    algorithm = new MOEADBuilder(problem, MOEADBuilder.Variant.MOEADMiriam)
            .setCrossover(crossover)
            .setMutation(mutation)
            .setMaxEvaluations(200000)
            .setPopulationSize(300)
            .setResultPopulationSize(300)
            .setNeighborhoodSelectionProbability(0.9)
            .setMaximumNumberOfReplacedSolutions(2)
            .setNeighborSize(20)
            .setFunctionType(AbstractMOEAD.FunctionType.TCHE)
            .build() ;
            
    JMetalLogger.logger.info("Run MultiPharmMRF");
    

    
    AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm)
        .execute() ;

    List<DoubleSolution> population = algorithm.getResult() ;
    long computingTime = algorithmRunner.getComputingTime() ;

    String[] temp = MultiPharmMRF.queryPath.split(".+?/(?=[^/]+$)");
    String queryBasename = (temp[temp.length-1].split("\\.(?=[^\\.]+$)"))[0];
    String[] temp2 = MultiPharmMRF.targetPath.split(".+?/(?=[^/]+$)");
    String targetBasename = (temp2[temp2.length-1].split("\\.(?=[^\\.]+$)"))[0];
    JMetalLogger.logger.info("Total execution time of "+targetBasename+": " + computingTime + "ms");
    System.out.println("Total execution time of "+targetBasename+": " + computingTime + "ms");

    //printFinalSolutionSet(population);
    //printFinalSolutionSet(population, algorithm.getName());
    String folderName = algorithm.getName()+"_"+queryBasename;
    printFinalSolutionSet(population, folderName, targetBasename);
    if (!referenceParetoFront.equals("")) {
      printQualityIndicators(population, referenceParetoFront) ;
    }
  }
}
