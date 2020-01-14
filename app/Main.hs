module Main where

import Lib
import Bio.Paf
import Bio.Paf.IO

main :: IO ()
main = do
    paf <- readPafFile "test/data/test.paf"
    ovlp <- pafToOvlp paf
    print $ take 5 ovlp
    (graph, nodeFromVertex, vertexFromKey) = graphFromEdges edgeList
    print nodeFromVertex
    print vertexFromkey
    print graph
  
