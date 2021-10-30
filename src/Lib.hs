{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE DeriveGeneric       #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Lib
    ( main
    ) where


import           Data.Text hiding 
 ( head, tail, zip, last, map, take
 )
import qualified Data.Text as T
import qualified Data.Text.IO as TIO
import           System.IO
import           GHC.Generics 
import qualified RNA
import qualified Data.List as L

basePairs :: Text
basePairs = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"

countBasePairs :: Text -> IO ()
countBasePairs bp = print $ go bp (0,0,0,0)
  where
   go :: Text -> (Int, Int, Int, Int) -> (Int, Int, Int, Int)
   go !"" !(a,b,c,d)      = (a,b,c,d)
   go !bps !(a,b,c,d)
    | (T.head bps) == 'A' = go (T.tail bps) (a+1, b, c, d)
    | (T.head bps) == 'G' = go (T.tail bps) (a, b+1, c, d)
    | (T.head bps) == 'T' = go (T.tail bps) (a, b, c+1, d)
    | (T.head bps) == 'C' = go (T.tail bps) (a, b, c, d+1)
    
    
dnaToRna :: Text -> Text
dnaToRna = T.replace "T" "U"

findSubStringPos :: Text -> Text -> [Int]
findSubStringPos base match = go base [] 0
  where
   matchLen = T.length match
   go t l i
    | T.length t >= matchLen = 
       case ((T.take matchLen t) == match) of
        True -> go (T.drop 1 t) (l ++ [i+1]) (i+1)
        False -> go (T.drop 1 t) (l) (i+1)
    | otherwise              = l

compliment :: Char -> Char
compliment 'A' = 'T' 
compliment 'T' = 'A' 
compliment 'C' = 'G' 
compliment 'G' = 'C' 
compliment _ = error "No Compliment"

complimentDna :: Text -> Text 
complimentDna = T.foldr (\x y -> T.cons (compliment x) y) "" 

rabbits :: Int -> Int -> Int
rabbits months litter = go 0 []
  where go m l
         | m == months = last l 
         | m == 0      = go (m+1) [1]
         | m == 1      = go (m+1) (l ++ [1])
         | otherwise   = go (m+1) (l ++ [((l !! (m-2))*litter) + (l !! (m-1))])

data FastaData = FastaData 
  { identifier :: Text 
  , bps        :: Text 
  } deriving (Eq, Show, Generic)

computeGc :: Fractional a => FastaData -> a
computeGc fd = (cmptd / fdLen)
  where fdLen = fromIntegral $ T.length (bps fd)
        cmptd = fromIntegral $ T.foldr (\x y -> if (x == 'G' || x == 'C') then (y + 1) else y) 0 (bps fd)

gcContent :: IO [FastaData]
gcContent = withFile "./gc.txt" ReadMode  (\h -> getLines h [] Nothing) 
  where 
    getLines h l fd = do
     e <- hIsEOF h
     case e of
      True -> 
       case fd of
        Just f -> return (l ++ [f])
        Nothing -> return l
      False -> do
       li <- TIO.hGetLine h
       case (T.head li) of
        '>' ->
         case (fd) of 
          Just fd -> getLines h (l ++ [fd]) (Just (FastaData (T.tail li) ""))
          Nothing -> getLines h l (Just (FastaData (T.tail li) ""))
        _ ->
         case (fd) of 
          Just fds -> getLines h l (Just (fds { bps = (append (bps fds) li) }))
          _ -> error "Invalid Format"

hammingDistance :: Text -> Text -> Int
hammingDistance t1 t2 = 
 Prelude.foldr (\(b1, b2) y -> if b1 /= b2 then (y+1) else y) 0 $ T.zip t1 t2


data Phenotype = HomoD
                | Hetero
                | HomoR
                deriving (Eq, Show, Generic)

inheritanceChances :: Fractional a => Phenotype -> Phenotype -> a
inheritanceChances HomoD _       = 1.0
inheritanceChances _ HomoD       = 1.0
inheritanceChances Hetero Hetero = 0.75
inheritanceChances HomoR Hetero  = 0.50
inheritanceChances Hetero HomoR  = 0.50
inheritanceChances HomoR HomoR   = 0.0

{-
probabilityOfInheritance :: Fractional a => Int -> Int -> Int -> a
probabilityOfInheritance hmd hz hmr = 
 ((fromIntegral hmd)/(fromIntegral tot)) + ((fromIntegral hz)/(fromIntegral tot))
 where tot = hmd + hz + hmr
-}
probabilityOfInheritance :: Fractional a => Int -> Int -> Int -> a
probabilityOfInheritance hmd hz hmr = (sum combs) / (fromIntegral $ Prelude.length combs)
 where 
  combs = map (\(x, y) -> inheritanceChances x y) $ combinations $ (take hmd (repeat HomoD)) ++ (take hz (repeat Hetero)) ++ (take hmr (repeat HomoR))

combinations :: [a] -> [(a,a)]
combinations vals = go vals []
  where 
   go :: [a] -> [(a,a)] -> [(a,a)]
   go [] l     = l
   go (a:as) l = go (as) (l ++ (Prelude.foldr (\x y -> y ++ [(a, x)]) [] as))

parseFastaData :: String -> IO [FastaData]
parseFastaData s = withFile s ReadMode  (\h -> getLines h [] Nothing) 
  where 
    getLines h l fd = do
     e <- hIsEOF h
     case e of
      True -> 
       case fd of
        Just f -> return (l ++ [f])
        Nothing -> return l
      False -> do
       li <- TIO.hGetLine h
       case (T.head li) of
        '>' ->
         case (fd) of 
          Just fd -> getLines h (l ++ [fd]) (Just (FastaData (T.tail li) ""))
          Nothing -> getLines h l (Just (FastaData (T.tail li) ""))
        _ ->
         case (fd) of 
          Just fds -> getLines h l (Just (fds { bps = (append (bps fds) li) }))
          _ -> error "Invalid Format"

updateMatrix :: Int -> Int -> a -> [[a]] -> [[a]]
updateMatrix row col val mat =  
  h ++ [(hc ++ [val] ++ (Prelude.tail tc))] ++ (Prelude.tail t)
  where (h, t) = Prelude.splitAt (row) mat
        (hc, tc) = Prelude.splitAt (col) (Prelude.head t)

getMatrix :: Int -> Int -> [[a]] -> a
getMatrix row col mat = ((!!) ((!!) mat row) col)

profileMatrix :: [Text] -> [[Int]]
profileMatrix t = go 0 (Prelude.replicate 4 (Prelude.replicate minLen 0))
  where 
   minLen = Prelude.minimum (map T.length t)
   arrLen = Prelude.length t
   go ind mat 
    | ind == minLen = mat
    | otherwise     = go (ind + 1) (addChar' 0 ind mat)

   addChar' r idx m
    | r >= arrLen = m
    | otherwise  = addChar' (r + 1) idx (addChar (T.head $ T.drop idx ((!!) t r)) m idx)

   addChar c m ind
    | c == 'A' = (updateMatrix 0 ind ((getMatrix 0 ind m) + 1) m)
    | c == 'C' = (updateMatrix 1 ind ((getMatrix 1 ind m) + 1) m)
    | c == 'G' = (updateMatrix 2 ind ((getMatrix 2 ind m) + 1) m)
    | c == 'T' = (updateMatrix 3 ind ((getMatrix 3 ind m) + 1) m)

printProfileMatrix :: [[Int]] -> IO ()
printProfileMatrix mat = do
  putStr "A: "
  mapM_ putStr $ L.intersperse (" "::String) (map show ((!!) mat 0))
  putStrLn ""
  putStr "C: "
  mapM_ putStr $ L.intersperse (" "::String) (map show ((!!) mat 1))
  putStrLn ""
  putStr "G: "
  mapM_ putStr $ L.intersperse (" "::String) (map show ((!!) mat 2))
  putStrLn ""
  putStr "T: "
  mapM_ putStr $ L.intersperse (" "::String) (map show ((!!) mat 3))
  putStrLn ""


matrixGetCol :: [[a]] -> Int -> [a]
matrixGetCol m i = map (\x -> Prelude.head $ Prelude.drop i x) m

maxIdx :: [Int] -> Int
maxIdx xs = Prelude.head $ Prelude.filter ((== Prelude.maximum xs) . (xs !!)) [0..] 

intToCodon :: Int -> Char
intToCodon i
 | i == 0 = 'A'
 | i == 1 = 'C'
 | i == 2 = 'G'
 | i == 3 = 'T'
 | otherwise  = error "No CODON"


getConsensus :: [[Int]] -> [Char]
getConsensus mat = maxs 
  where a = (!!) mat 0
        c = (!!) mat 1
        g = (!!) mat 2
        t = (!!) mat 3
        minLen = Prelude.minimum (map Prelude.length mat)
        cols = map (\x -> matrixGetCol mat x) [0..(minLen-1)]
        maxs = map intToCodon (map maxIdx cols)
        {-
		colMax currMax idx
		 | idx == 4 = currMax
		 | otherwise = 
		  case (!!) mat idx
		-}

consProf :: IO ()
consProf = do
 p <- parseFastaData "./cons.txt"
 let s = map bps p
     pm = profileMatrix s
 printProfileMatrix pm
 print $ getConsensus pm
 print "Done"

mortalRabbits :: Int -> Int -> Int
mortalRabbits months lifespan = go 0 []
  -- where go m l
  where 
        go :: Int -> [Int] -> Int
        go m l
         | m == months = last l 
         | m == 0      = go (m+1) [1]
         | m == 1      = go (m+1) (l ++ [1])
         | otherwise   = go (m+1) (calcRabbits l m)
        
        calcRabbits :: [Int] -> Int -> [Int]
        calcRabbits e monthselapsed
         | monthselapsed <= lifespan = (e ++ [((e !! (monthselapsed-2))) + (e !! (monthselapsed-1))])
         | otherwise = (e ++ [((e !! (monthselapsed-2))) + (e !! (monthselapsed-1)) - (e !! (lifespan-1))])


overlapGraph :: [FastaData] -> Int -> [(Text, Text)]
overlapGraph fdl matchLen = Prelude.foldr (\x y -> matchAgainst x fdl y) [] fdl
  where 
   matchAgainst :: FastaData -> [FastaData] -> [(Text, Text)] -> [(Text, Text)]
   matchAgainst fd (x:xs) col
    | xs == [] = col 
    | (identifier fd) == (identifier x) = matchAgainst fd xs col
    | suf == (T.take matchLen (bps x)) = matchAgainst fd xs ([(identifier fd,identifier x)] ++ col)
    | otherwise = matchAgainst fd xs col

     where 
      suf = T.reverse $ T.take matchLen $ T.reverse $ bps fd
  

expectedOffspring :: Int -> Int -> Int -> Int -> Int -> Int -> Double
expectedOffspring a b c d e f = 2 * sum
 [ (fromIntegral a) * 1.00 
 , (fromIntegral b) * 1.00 
 , (fromIntegral c) * 1.00 
 , (fromIntegral d) * 0.75 
 , (fromIntegral e) * 0.25 
 , (fromIntegral f) * 0.00 
 ]

-- longestCommonSubstring :: [FastaData] -> [Text]
-- longestCommonSubstring fd =
--   where
--    go :: [Text] -> [Text] -> [Text]
--    go (x:xs) strs 
{-
matchSubStrs :: Text -> Text -> [Text]
matchSubStrs m1 m2 =
 where 
  go :: Text -> Text -> [Text] -> [Text]
  go "" m col = col
  go s fs 
-}

main :: IO ()
main = do
  print $ "RUNNING"
  {-
  countBasePairs basePairs
  print $ dnaToRna basePairs
  print $ complimentDna basePairs
  print $ rabbits 6 3
  gc <- gcContent
  print $ Prelude.foldr (\x (gc, n) -> if (computeGc x) > gc then ((computeGc x), identifier x) else (gc, n)) (0.0, "") gc
  print $ hammingDistance "GAGCCTACTAACGGGAT" "CATCGTAATGACGGCCT"
  print $ probabilityOfInheritance 2 2 2
  print $ map RNA.codonToChar $ RNA.parseRNAText "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
  print $ findSubStringPos "GATATATGCATATACTT" "ATAT"
  -}
  consProf
  print $ mortalRabbits 6 3
  graphData <- parseFastaData "./graph.txt"
  print $ overlapGraph graphData 3
  print $ expectedOffspring 1 0 0 1 0 1
