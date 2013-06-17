import sys
import partition as pa

test = pa.genotype(sys.argv[1])
test.read_file()
test.optiUberScan()
test.optiRLScan()
test.optiLRScan()
test.getoptiCore()
test.pesiUberScan()
test.pesiRLScan()
test.pesiLRScan()
test.getpesiCore()

print "The optimal Core is:"
test.displayoptiCore()
print "The pesitic Core is:"
test.displaypesiCore()

test.findoptiSusSNP()
print "optimal susSNP"
print test.opticriticalSNP
test.findpesiSusSNP()
print "pesimal susSNP"
print test.pesicriticalSNP


print "display opti CUber"
test.displayoptiCUber()
print "display pesi CUber"
test.displaypesiCUber()

test.assignoptiKGroup()
test.maximaloptiKCover()

print "opti maximal-k-cover"
test.displayoptiPath()

test.assignpesiKGroup()
test.maximalpesiKCover()

print "pesi maximal-k-cover"
test.displaypesiPath()
