{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}
module Lib
    ( pafToOvlp,
      ctgPairToEdge
    ) where

import Bio.Paf
import Bio.Paf.IO
import Data.Char
import Data.Graph
import Data.List
import qualified Data.Text as T
import Data.Text.Read 
import qualified Data.Vector as V
import GHC.Generics
import GHC.Int

data Contig = Contig {
   contigId :: Int,
   posStrand :: Bool
} deriving (Show)

contigAsInt :: Contig -> Int
contigAsInt a = case (posStrand a) of
  True -> (contigId a) * 2
  False -> (contigId a) * 2 + 1 

hangLength :: Int32
hangLength = 1000

adjacentLength :: Int32
adjacentLength = 10000

minMapQ :: Int32
minMapQ = 0

compareByQName :: Aln -> Aln -> Bool
compareByQName a b = (_qname a) == (_qname b)

isHangingAlignment:: Aln -> Bool
isHangingAlignment aln = case (_strand aln) of
    '+' -> (_tlen aln - hangLength) <= (_tend aln)
    '-' -> (_tstart aln) <= hangLength

isOverlapAlignment :: Aln -> Aln -> Bool
isOverlapAlignment a b = isHangingAlignment a && isHangingAlignment b && (_qstart b - _qend a) <= adjacentLength

extractOverlapAlignment :: [Aln] -> [(Contig, Contig)]
extractOverlapAlignment aln = map (uncurry convertToEdge) $ filter (uncurry isOverlapAlignment) $ map (\x -> (x !! 0, x !! 1)) $ eachCons 2 aln

pafToOvlp :: Paf -> [(Contig, Contig)]
pafToOvlp paf = do
    let grouped = groupBy compareByQName $ filter (\x -> _mapq x > minMapQ) $ V.toList $ _alns paf 
    let edgeList = concat $ map extractOverlapAlignment $ filter (\x -> (length x) >= 2) grouped 
    edgeList
    --let (graph, nodeFromVertex, vertexFromKey) = graphFromEdges edgeList
    --return graph

ctgPairToEdge :: (Contig, Contig) -> Edge
ctgPairToEdge (a, b) = (contigAsInt a, contigAsInt b) 

convertToEdge :: Aln -> Aln -> (Contig, Contig)
convertToEdge a b = (convertToContig a, convertToContig b)

convertToContig :: Aln -> Contig
convertToContig a = case (_strand a) of
    '+' -> Contig (contigIdToInt $ _tname a) True
    '-' -> Contig (contigIdToInt $ _tname a) False

contigIdToInt :: T.Text -> Int
contigIdToInt x = (myRead x) !! 0 

-- https://stackoverflow.com/questions/28506961/how-to-parse-text-and-extract-integer
myRead :: T.Text -> [Int]
myRead = map num                 -- Convert the remaining elements into Ints
       . filter (not . T.null)   -- Drop all empty words
       . map (T.filter isDigit)  -- Drop all non-digits in each word (including signs!)
       . T.words                 -- Chop the string into words

num :: T.Text -> Int 
num = either (error . show) fst  -- Throw an exception if it wasn't a signed decimal
    . signed decimal             -- Read a signed decimal

--Author: Mark D. Blackwell
--Change dates:
--(mdb) September 23, 2011 - create
--(mdb) September 27, 2011 - add zero and negative to eachCons 
eachCons :: Int -> [a] -> [[a]]
eachCons n x
    | n  < 0                   = error ("eachCons " ++ showsPrec 11 n "")
    | n  > (length $ take n x) = error ("eachCons " ++ showsPrec 11 n "") -- don't know how to put x in there.
    | n == 0 = []: [ [] | _ <- x] -- Add one, to fit the pattern.
    | n  > 0 = consecutives
    where
      sequences = map (`drop` x) [0..(n-1)]
      consecutives = [y | y <- transpose sequences, n==length y]
 
