module Lib
    ( pafToOvlp
    ) where

import Bio.Paf
import Bio.Paf.IO
import Data.Graph
import qualified Data.Text as T
import GHC.Generics

data Contig = Contig {
   contigId :: Int,
   strand :: Bool
} deriving (Generic, Show)

someFunc :: IO ()
someFunc = putStrLn "someFunc"

hangLength :: Int
hangLength = 1000

adjacentLength :: Int
adjacentLength = 10000

compareByQName :: Aln -> Aln -> Bool
compareByQName a b = (_qName a) == (_qname b)

isHangingAlignment:: Aln -> Bool
isHangingAlignment aln = case (_strand aln) of
    '+' -> (_tlen aln - hangLength) <= (_tend aln)
    '-' -> (_tstart) <= hangLength

isOverlapAlignment :: Aln -> Aln -> Bool
isOverlapAlignment a b = isHangingAlignment a && isHangingAlignment b && (_qstart b - _qend a) <= adjacentLength

extractOverlapAlignment :: [Aln] -> [(Contig, Contig)]
extractOverlapAlignment aln = map convertToEdge $ filter isOverlapAlignment $ eachCons aln

pafToOvlp :: Paf -> [Contig]
pafToOvlp paf = do
    let grouped = groupBy compareByQName $ filter (\x -> _mapq x > 0) paf 
    let edgeList = map (extractOverlapAlignment) grouped 
    return edgeList
    --let (graph, nodeFromVertex, vertexFromKey) = graphFromEdges edgeList
    --return graph

convertToEdge :: Aln -> Aln -> (Contig, Contig)
convertToEdge a b = (convetToContig a, convertToContig b)

convertToContig :: Aln -> Contig
convertToContig a = case (_strand a) of
    '+' -> Contig (contigIdToInt $ _tname a) True
    '-' -> Contig (contigIdToInt $_tname a) False

contigIdToInt :: T.Text -> int
contigIdToInt x = (myRead x) !! 0 

myRead :: T.Text -> [Int]
myRead = map num                 -- Convert the remaining elements into Ints
       . filter (not . T.null)   -- Drop all empty words
       . map (T.filter isDigit)  -- Drop all non-digits in each word (including signs!)
       . T.words                 -- Chop the string into words

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
 
