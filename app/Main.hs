module Main where

import Lib
import Bio.Paf
import Bio.Paf.IO
import qualified Data.Vector as V
import Data.List
import System.Environment (getArgs)
import Data.Graph

compareByQName :: Aln -> Aln -> Bool
compareByQName a b = (_qname a) == (_qname b)

fileName :: [String] -> String
fileName args = do
    case 1 <= length args of
        True -> args !! 0
        False -> "test/data/test.paf"

main :: IO ()
main = do
    args <- getArgs
    let paf = fileName args
    
    paf <- readPafFile paf
--    let grouped = filter (\x -> (length x) >= 2) $groupBy compareByQName $ filter (\x -> _mapq x > 0) $ V.toList $ _alns paf 
--    print $ grouped
    let ovlp = pafToOvlp paf
    mapM_ print ovlp
    --print $ ovlp
    let edges = map ctgPairToEdge ovlp
    let maxEdgeId = maximum $ concat $ map pairToList edges
    let blackEdges = map (\x -> (x!!0, x!!1)) $ slices 2 $ [0..maxEdgeId] 
    let totalEdges = edges `mappend` blackEdges
    print $ totalEdges
    let graph = buildG (0, maxEdgeId) $ totalEdges
    let start = Contig 10 False
    let stop = Contig 15 True
    print $ path graph (contigAsInt start) (contigAsInt stop)
    --print $ dfs graph [contigAsInt start]
    --let start2 = Contig 10 True
    --print $ dfs graph [contigAsInt start2]
    --print graph
    --let (graph, nodeFromVertex, vertexFromKey) = graphFromEdges ovlp
    --print nodeFromVertex
    --print vertexFromKey
    --print graph

slices n s | length s < n = [] | otherwise = (take n s):slices n (drop n s)

pairToList :: (a, a) -> [a]
pairToList (x,y) = [x,y]
