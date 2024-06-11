# MultiPharm

Source code of MultiPharm under Mozilla Public License Version 2.0, https://www.mozilla.org/en-US/MPL/2.0/
For non-academics: a licence needed. Contact with the authors.

jMetal and Openeye Libraries belong to their authors and access to them depends on any applicable restrictions.
Additionally, a valid OpenEye license is mandatory to run the software.

For more information about libraries go to https://github.com/jMetal/jMetal and https://www.eyesopen.com/.

# How to run
1. mvn package
2. export CLASSPATH=jmetal-core/target/jmetal-core-5.8-jar-with-dependencies.jar:jmetal-problem/target/jmetal-problem-5.8-jar-with-dependencies.jar:jmetal-exec/target/jmetal-exec-5.8-jar-with-dependencies.jar:jmetal-problem/target/jmetal-problem-5.8-jar-with-dependencies.jar
3. java org.uma.jmetal.runner.multiobjective.MOEADRunnerMultiPharm org.uma.jmetal.problem.multiobjective.MultiPharm [path-to-query] [path-to-target]
4. Results are saved in a folder with the name of the query molecule.
