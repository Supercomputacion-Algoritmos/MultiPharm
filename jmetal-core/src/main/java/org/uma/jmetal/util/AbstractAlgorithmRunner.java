package org.uma.jmetal.util;

import org.uma.jmetal.qualityindicator.impl.*;
import org.uma.jmetal.qualityindicator.impl.hypervolume.PISAHypervolume;
import org.uma.jmetal.solution.Solution;
import org.uma.jmetal.util.fileoutput.SolutionListOutput;
import org.uma.jmetal.util.fileoutput.impl.DefaultFileOutputContext;
import org.uma.jmetal.util.front.Front;
import org.uma.jmetal.util.front.imp.ArrayFront;
import org.uma.jmetal.util.front.util.FrontNormalizer;
import org.uma.jmetal.util.front.util.FrontUtils;
import org.uma.jmetal.util.point.PointSolution;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

/**
 * Abstract class for Runner classes
 *
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 */
public abstract class AbstractAlgorithmRunner {
  /**
   * Write the population into two files and prints some data on screen
   * @param population
   */
  public static void printFinalSolutionSet(List<? extends Solution<?>> population) {

    new SolutionListOutput(population)
        .setSeparator("\t")
        .setVarFileOutputContext(new DefaultFileOutputContext("VAR.tsv"))
        .setFunFileOutputContext(new DefaultFileOutputContext("FUN.tsv"))
        .print();

    JMetalLogger.logger.info("Random seed: " + JMetalRandom.getInstance().getSeed());
    JMetalLogger.logger.info("Objectives values have been written to file FUN.tsv");
    JMetalLogger.logger.info("Variables values have been written to file VAR.tsv");
  }
  
  public static void printFinalSolutionSet(List<? extends Solution<?>> population, String namef) {

    new SolutionListOutput(population)
        .setSeparator("\t")
        .setVarFileOutputContext(new DefaultFileOutputContext(namef+"/VAR.tsv"))
        .setFunFileOutputContext(new DefaultFileOutputContext(namef+"/FUN.tsv"))
        .print();

    JMetalLogger.logger.info("Random seed: " + JMetalRandom.getInstance().getSeed());
    JMetalLogger.logger.info("Objectives values have been written to file FUN.tsv");
    JMetalLogger.logger.info("Variables values have been written to file VAR.tsv");
  }
  
  public static void printFinalSolutionSet(List<? extends Solution<?>> population, String namef, String ending) {

    // Create the folder if it doesnt exist
    boolean successful = false;
    int counter = 1;
    String newFolder = namef +"_" +counter;
    File dir = new File(newFolder);
    while (!successful){
      successful = dir.mkdir();
      if (successful)
      {
        // creating the directory succeeded
        System.out.println("directory "+newFolder+ " was created successfully");
      }
      else
      {
        counter += 1;
        newFolder = namef +"_" +counter;
        dir = new File(newFolder);
      }
    }
    namef = newFolder;
    new SolutionListOutput(population)
        .setSeparator("\t")
        .setVarFileOutputContext(new DefaultFileOutputContext(namef+"/VAR_"+ending+".tsv"))
        .setFunFileOutputContext(new DefaultFileOutputContext(namef+"/FUN_"+ending+".tsv"))
        .print();

    //JMetalLogger.logger.info("Random seed: " + JMetalRandom.getInstance().getSeed());
    //JMetalLogger.logger.info("Objectives values have been written to file FUN.tsv");
    //JMetalLogger.logger.info("Variables values have been written to file VAR.tsv");
  }

  /**
   * Print all the available quality indicators
   * @param population
   * @param paretoFrontFile
   * @throws FileNotFoundException
   */
  public static <S extends Solution<?>> void printQualityIndicators(List<S> population, String paretoFrontFile)
      throws FileNotFoundException {
    Front referenceFront = new ArrayFront(paretoFrontFile);
    FrontNormalizer frontNormalizer = new FrontNormalizer(referenceFront) ;

    Front normalizedReferenceFront = frontNormalizer.normalize(referenceFront) ;
    Front normalizedFront = frontNormalizer.normalize(new ArrayFront(population)) ;
    List<PointSolution> normalizedPopulation = FrontUtils
        .convertFrontToSolutionList(normalizedFront) ;

    String outputString = "\n" ;
    outputString += "Hypervolume (N) : " +
        new PISAHypervolume<PointSolution>(normalizedReferenceFront).evaluate(normalizedPopulation) + "\n";
    outputString += "Hypervolume     : " +
        new PISAHypervolume<S>(referenceFront).evaluate(population) + "\n";
    outputString += "Epsilon (N)     : " +
        new Epsilon<PointSolution>(normalizedReferenceFront).evaluate(normalizedPopulation) +
        "\n" ;
    outputString += "Epsilon         : " +
        new Epsilon<S>(referenceFront).evaluate(population) + "\n" ;
    outputString += "GD (N)          : " +
        new GenerationalDistance<PointSolution>(normalizedReferenceFront).evaluate(normalizedPopulation) + "\n";
    outputString += "GD              : " +
        new GenerationalDistance<S>(referenceFront).evaluate(population) + "\n";
    outputString += "IGD (N)         : " +
        new InvertedGenerationalDistance<PointSolution>(normalizedReferenceFront).evaluate(normalizedPopulation) + "\n";
    outputString +="IGD             : " +
        new InvertedGenerationalDistance<S>(referenceFront).evaluate(population) + "\n";
    outputString += "IGD+ (N)        : " +
        new InvertedGenerationalDistancePlus<PointSolution>(normalizedReferenceFront).evaluate(normalizedPopulation) + "\n";
    outputString += "IGD+            : " +
        new InvertedGenerationalDistancePlus<S>(referenceFront).evaluate(population) + "\n";
    outputString += "Spread (N)      : " +
        new Spread<PointSolution>(normalizedReferenceFront).evaluate(normalizedPopulation) + "\n";
    outputString += "Spread          : " +
        new Spread<S>(referenceFront).evaluate(population) + "\n";
//    outputString += "R2 (N)          : " +
//        new R2<List<DoubleSolution>>(normalizedReferenceFront).runAlgorithm(normalizedPopulation) + "\n";
//    outputString += "R2              : " +
//        new R2<List<? extends Solution<?>>>(referenceFront).runAlgorithm(population) + "\n";
    outputString += "Error ratio     : " +
        new ErrorRatio<List<? extends Solution<?>>>(referenceFront).evaluate(population) + "\n";
    
    JMetalLogger.logger.info(outputString);
  }
}
