import sys
import partition as pa

test = pa.haplotype(sys.argv[1])
test.read_file()
test.UberScan()
test.RLScan()
test.LRScan()
test.getCore()

print "displayCore"
test.displayCore()

test.findSusSNP()
print "susSNP"
print test.criticalSNP


print "displayUber"
test.displayCUber()

test.assignKGroup()
test.maximalKCover()

print "maximal-k-cover"
test.displayPath()




