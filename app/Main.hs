module Main where

import Lib
import Bio.Paf
import Bio.Paf.IO
import qualified Data.Vector as V
import Data.List
import Data.Graph

compareByQName :: Aln -> Aln -> Bool
compareByQName a b = (_qname a) == (_qname b)


main :: IO ()
main = do
    paf <- readPafFile "test/data/test.paf"
    let grouped = filter (\x -> (length x) >= 2) $groupBy compareByQName $ filter (\x -> _mapq x > 0) $ V.toList $ _alns paf 
    print $ grouped
    let ovlp = pafToOvlp paf
    print $ take 5 ovlp
    let edges = map ctgPairToEdge ovlp
    let graph = buildG (0, maximum $ concat $ map pairToList edges) $ edges
    print graph
    --let (graph, nodeFromVertex, vertexFromKey) = graphFromEdges ovlp
    --print nodeFromVertex
    --print vertexFromKey
    --print graph

pairToList :: (a, a) -> [a]
pairToList (x,y) = [x,y]
