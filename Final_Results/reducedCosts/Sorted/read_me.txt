Performance Profile by D. Salvagnin (2016)
Internal use only, not to be distributed   

Dolan ED, Moŕe J (2002) Benchmarking optimization software with performance profiles. Mathematical Programming 91(2):201–213

@article{dolan2002benchmarking,
  title={Benchmarking optimization software with performance profiles},
  author={Dolan, E. D. and Mor{\'e}, J.J.},
  journal={Mathematical Programming},
  volume={91},
  number={2},
  pages={201--213},
  year={2002},
  publisher={Springer}
}


Usage:

python ../perfprof.py -D , -T 3600 -S 2 -M 20 lagr.csv pp.pdf -P "all instances, shift 2 sec.s"  


Parameters:

self.parser.add_option("-D", "--delimiter", dest="delimiter", default=None, help="delimiter for input files")
self.parser.add_option("-M", "--maxratio", dest="maxratio", default=4, type=int, help="maxratio for perf. profile")
self.parser.add_option("-S", "--shift", dest="shift", default=0, type=float, help="shift for data")
self.parser.add_option("-L", "--logplot", dest="logplot", action="store_true", default=False, help="log scale for x")
self.parser.add_option("-T", "--timelimit", dest="timelimit", default=1e99, type=float, help="time limit for runs")
self.parser.add_option("-P", "--plot-title", dest="plottitle", default=None, help="plot title")
self.parser.add_option("-X", "--x-label", dest="xlabel", default='Time Ratio', help="x axis label")
self.parser.add_option("-B", "--bw", dest="bw", action="store_true", default=False, help="plot B/W")



example:
python3 ./perfprof.py -D , -T 3600 -S 2 -M 10 GenericConcorde_VS_Generic_VS_Legacy.csv  GenericConcorde_VS_Generic_VS_Legacy.pdf -P "all instances, shift 2 sec.s"


python3 ./perfprof.py -D , -T 3600 -S 2 -M 10 Loop_gap50_20_10_5_0.csv  Loop_gap50_20_10_5_0.pdf -P "all instances, shift 2 sec.s"

python3 ./perfprof.py -D , -T 3600 -S 2 -M 10 Concorde.csv  Concorde.pdf -P "all instances, shift 2 sec.s"




python3 ./perfprof.py -D , -T 3600 -S 2 -M 10 newCompact.csv Concorde.pdf -P "Compact models"

python3 ./perfprof.py -D , -T 600 -S 2 -M 12 newLoop.csv newLoop.pdf -P "Loop variants"

python3 ./perfprof.py -D , -T 600 -S 2 -M 95 newCompactVsLoop.csv newCompactVsLoop.pdf -P "Compact vs Loop"

python3 ./perfprof.py -D , -T 600 -S 2 -M 10 newLegacyVsGeneric.csv newLegacyVsGeneric.pdf -P "Legacy vs Generic"

python3 ./perfprof.py -D , -T 600 -S 2 -M 10 newLoopVsCallback.csv newLoopVsCallback.pdf -P "Loop vs Callback"

python3 ./perfprof.py -D , -T 600 -S 2 -M 1.5 newConcorde.csv newConcorde.pdf -P "Concorde versions"

python3 ./perfprof.py -D , -T 600 -S 2 -M 1.10 newConcordeVsGeneric.csv newConcordeVsGeneric.pdf -P "Concorde"


python3 ./perfprof.py -D , -M 1.1 newConstructiveHeuristics.csv newConstructiveHeuristics.pdf -P "Constructive Heuristics"

python3 ./perfprof.py -D , -M 2 newConstructiveHeurisricsLong.csv newConstructiveHeurisricsLong.pdf -P "Constructive Heuristics"

python3 ./perfprof.py -D , -M 1.05 newRefiningHeuristics.csv newRefiningHeuristics.pdf -P "Standalone Refining Heuristics"

python3 ./perfprof.py -D , -M 1.1 newMatheuristics.csv newMatheuristics.pdf -P "Matheuristics"

python3 ./perfprof.py -D , -M 1.1 newAllHeuristics.csv newAllHeuristics.pdf -P "Heuristics"

python3 ./perfprof.py -D , -M 1.1 newHeuristisNo2Opt.csv newHeuristisNo2Opt.pdf -P "Constructive Heuristics"

python3 ./perfprof.py -D , -M 1.1 newHeuristisNo2Opt.csv newHeuristisNo2Opt.pdf -P "Constructive Heuristics"

python3 ./perfprof.py -D , -M 1.1 oldVsNewPatching.csv oldVsNewPatching.pdf -P "oldVsNewPatching"

python3 ./perfprof.py -D , -M 1.3 new2Opt.csv new2Opt.pdf -P "Constructive Heuristics with or without 2-opt"

python3 ./perfprof.py -D , -M 1.1 newCPLEXBasedHeuristics.csv newCPLEXBasedHeuristics.pdf -P "CPLEX-based Heuristics"

python3 ./perfprof.py -D , -M 1.1 newFinalComparison.csv newFinalComparison.pdf -P "Heuristics"

python3 ./perfprof.py -D , -M 1.1 newPatchingPolishing.csv newPatchingPolishing.pdf -P "Patching + RINS and Polishing"

python3 ./perfprof.py -D , -M 1.1 newExactSolution.csv newExactSolution.pdf -P "Comparison with CPLEX generic callback and its heuristics"

python3 ./perfprof.py -D , -T 3600 -S 2 -M 10 1.csv 1.pdf -P "Time comparison par = 1.5 lookahead = 50" -X "Time Ratio"

python3 ./perfprof.py -D , -S 100 -M 10 2.csv 2.pdf -P "Node Comparison par = 1.5 lookahead = 50" -X "Node Ratio"
