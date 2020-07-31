{-# LANGUAGE DeriveGeneric, DeriveAnyClass, OverloadedStrings #-}
module Main where

import Lib
import Bio.Paf
import Bio.Paf.IO
import Data.Aeson.Text (encodeToLazyText)
import Data.Aeson (ToJSON)
import qualified Data.Vector as V
import Data.List
import System.Environment (getArgs)
import Data.Graph
import Data.Graph.Inductive hiding (dfs, Graph)
import Data.Graph.Inductive.Example (genUNodes)
import Data.Maybe (fromMaybe)
import Data.Text.Lazy.IO as I

--import Data.Graph.Inductive.Graph
--import Data.Graph.Inductive.Tree

compareByQName :: Aln -> Aln -> Bool
compareByQName a b = (_qname a) == (_qname b)

fileName :: [String] -> String
fileName args = do
    case 1 <= length args of
        True -> args !! 0
        False -> "test/data/test.paf"

countFreqEdges :: [(Node, Node)] -> [(Node, Node, Int)]
countFreqEdges = elemCountOnce 

--https://qiita.com/mather314/items/dbae21e2fcde2ea74669
elemCountOnceOld :: Eq a => [a] -> [(a, Int)]
elemCountOnceOld [] = []
elemCountOnceOld (x:xs) = (x, (length $ filter (== x) xs) + 1) : elemCountOnceOld (filter (/= x) xs)

elemCountOnce :: Eq a => [(a, a)] -> [(a, a, Int)]
elemCountOnce [] = []
elemCountOnce ((a,b):xs) = (a, b, (length $ filter (== (a, b)) xs) + 1) : elemCountOnce (filter (/= (a, b)) xs)

createGraph :: [(Contig, Contig)] -> Graph
createGraph ovlp = do
    let edges = concat $ map ctgPairToEdges ovlp
    let maxEdgeId = maximum $ concat $ map pairToList edges
    buildG (0, maxEdgeId) $ edges

createInductiveGraph :: [(Contig, Contig)] -> Gr Contig Int -- Edge label is frequency, and node label is currently none
createInductiveGraph ovlp = do
    let edges = concat $ map ctgPairToEdges ovlp
    let maxEdgeId = maximum $ concat $ map pairToList edges
    --mkGraph (genUNodes maxEdgeId) $ countFreqEdges edges 
    mkGraph (genNodes maxEdgeId) $ countFreqEdges edges 

edgeDump :: [(Contig, Contig)] -> String -> IO ()
edgeDump ovlp json = do
    --let edges = concat $ map ctgPairToEdges ovlp
    I.writeFile json (encodeToLazyText $ elemCountOnceOld ovlp)

edgeDumpJson :: [(Contig, Contig)] -> [String] -> IO ()
edgeDumpJson ovlp args = do
    case 2 <= length args of
        True -> edgeDump ovlp (args !! 1)
        False -> return ()

main :: IO ()
main = do
    args <- getArgs
    let paf = fileName args
    
    --paf <- readPafFile paf
--    let grouped = filter (\x -> (length x) >= 2) $groupBy compareByQName $ filter (\x -> _mapq x > 0) $ V.toList $ _alns paf 
--    print $ grouped
    --let ovlp = pafToOvlp paf
    ovlp <- pafToOvlp <$> readPafFile paf
    edgeDumpJson ovlp args
    --mapM_ print ovlp
    --print $ ovlp
    --let blackEdges = map (\x -> (x!!0, x!!1)) $ slices 2 $ [0..maxEdgeId] 
    --let totalEdges = edges `mappend` blackEdges
    --print $ totalEdges
    --let edges = concat $ map ctgPairToEdges ovlp
    --let maxEdgeId = maximum $ concat $ map pairToList edges
    let graph = createGraph ovlp --buildG (0, maxEdgeId) $ edges
    let start = Contig 15 "" True -- Contig 15 is forward strand at the left side of chr19.
    let stop = Contig 10 "" True -- Contig 10 is forward strand at the right side of chr19, so the end strand should be forward.
    let start_ = contigAsInt start
    let stop_ = contigAsInt stop 
    print $ path graph (contigAsInt start) (contigAsInt stop)
    print $ dfs graph [contigAsInt start]
    let igraph = createInductiveGraph ovlp
    print "Onehop from start/stop"
    print $ context igraph start_ --igraph
    print $ lsuc igraph start_ --igraph
    print $ out igraph start_ --igraph
    print $ lpre igraph stop_ --igraph
    print $ inn igraph stop_ --igraph
    print "Shortestpath"
    print $ sp (contigAsInt start) (contigAsInt stop) igraph
    print "Other statistics"
    print $ lesp (contigAsInt start) (contigAsInt stop) igraph -- the minimum
    print $ spLength (contigAsInt start) (contigAsInt stop) igraph -- length of nodes
    prettyPrint $ subgraph [(contigAsInt start), (contigAsInt stop)] igraph
    print "Filtered shortestpath"
    let dygraph = efilter (\x -> getThd3 x >= 2) igraph  -- Filter out the edges greater than 
    print $ sp (contigAsInt start) (contigAsInt stop) dygraph
    print $ lesp (contigAsInt start) (contigAsInt stop) dygraph
    print "Merged Voronoiset"
    let fromStart = voronoiSet start_ $ gvdOut [start_] igraph
    let toStart = voronoiSet stop_ $ gvdIn [stop_] igraph
    print $ intersect fromStart toStart
    --prettyPrint $ subgraph (intersect fromStart toStart) igraph
    let path = fromMaybe [] $ sp (contigAsInt start) (contigAsInt stop) igraph
    print "Path Degree"
    print $ map (outdeg igraph) path
    print $ map (indeg igraph) path
    print "Weight-inverted Graph"
    let dyngraph = emap (\e -> (1.0 / fromIntegral e)) igraph
    print $ sp (contigAsInt start) (contigAsInt stop) dyngraph
    print $ lesp (contigAsInt start) (contigAsInt stop) dyngraph

    let chr1_1 = Contig 466 "" False -- Contig 15 is forward strand at the left side of chr19.
    let chr1_2 = Contig 112 "" False -- Contig 15 is forward strand at the left side of chr19.
    let chr5_1 = Contig 111 "" False -- Contig 15 is forward strand at the left side of chr19.

    print "edge for unique contig"
    print $ sp (contigAsInt start) (contigAsInt chr1_1) dyngraph
    print $ sp (contigAsInt start) (contigAsInt chr1_2) dyngraph
    print $ sp (contigAsInt start) (contigAsInt chr5_1) dyngraph
    print $ sp (contigAsInt chr5_1) (contigAsInt stop) dyngraph
    print $ sp (contigAsInt chr1_1) (contigAsInt stop) dyngraph
    print $ sp (contigAsInt chr1_2) (contigAsInt stop) dyngraph

    print "edge for unique contig w/context"
    print $ map (context igraph) $ fromMaybe [] $ sp (contigAsInt start) (contigAsInt stop) dyngraph
    print $ map (context igraph) $ fromMaybe [] $ sp (contigAsInt start) (contigAsInt chr1_1) dyngraph
    
    print "Onehop from start/stop"
    print $ context igraph 936 --igraph
    print $ map (context igraph) $ neighbors' $ context igraph 936
    print $ map (context igraph) $ neighbors' $ context igraph (contigAsInt chr1_1)
    print $ map (context igraph) $ neighbors' $ context igraph (contigAsInt chr1_2)
    print $ map (context igraph) $ neighbors' $ context igraph (contigAsInt chr5_1)
    --let start2 = Contig 10 True
    --print $ dfs graph [contigAsInt start2]
    --print graph
    --let (graph, nodeFromVertex, vertexFromKey) = graphFromEdges ovlp
    --print nodeFromVertex
    --print vertexFromKey
    --print graph

getThd3 (_,_,c)=c
slices n s | length s < n = [] | otherwise = (take n s):slices n (drop n s)

pairToList :: (a, a) -> [a]
pairToList (x,y) = [x,y]
