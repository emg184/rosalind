{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}

module RNA where


import           Data.Text hiding 
 ( head, tail, zip, last, map, take
 )
import qualified Data.Text as T
import qualified Data.Text.IO as TIO
import           System.IO
import           GHC.Generics 

parseRNAText :: Text -> [RNACodon]
parseRNAText t = go t []
  where
   go bps l
    | T.length bps >= 3 = go (T.drop 3 bps) (l ++ [tripletToRNACodon $ T.take 3 bps])
    | otherwise         = l



data RNACodon =
   Ala
 | Arg
 | Asn
 | Asp
 | Gly
 | Ser
 | Trp
 | Opal
 | Cys
 | Tyr
 | Ochre
 | Amber
 | His
 | Gln
 | Lys
 | Glu
 | Pro
 | Thr
 | Val
 | Met
 | Ile
 | Leu
 | Phe
 deriving (Eq, Show, Generic)

tripletToRNACodon :: Text -> RNACodon
tripletToRNACodon t
  | t == "UUU" = Phe
  | t == "UUC" = Phe
  | t == "UUA" = Leu
  | t == "UUG" = Leu
  | t == "CUU" = Leu
  | t == "CUC" = Leu
  | t == "CUA" = Leu
  | t == "CUG" = Leu
  | t == "AUU" = Ile
  | t == "AUC" = Ile
  | t == "AUA" = Ile
  | t == "AUG" = Met
  | t == "GUU" = Val
  | t == "GUC" = Val
  | t == "GUA" = Val
  | t == "GUG" = Val
  | t == "GUG" = Val
  | t == "UCU" = Ser
  | t == "UCC" = Ser
  | t == "UCA" = Ser
  | t == "UCG" = Ser
  | t == "CCU" = Pro
  | t == "CCC" = Pro
  | t == "CCA" = Pro
  | t == "CCG" = Pro
  | t == "ACU" = Thr
  | t == "ACC" = Thr
  | t == "ACA" = Thr
  | t == "ACG" = Thr
  | t == "GCU" = Ala
  | t == "GCC" = Ala
  | t == "GCA" = Ala
  | t == "GCG" = Ala
  | t == "UAU" = Tyr
  | t == "UAC" = Tyr
  | t == "UAA" = Ochre
  | t == "UAG" = Amber
  | t == "CAU" = His
  | t == "CAC" = His
  | t == "CAA" = Gln
  | t == "CAG" = Gln
  | t == "AAU" = Asn
  | t == "AAC" = Asn
  | t == "AAA" = Lys
  | t == "AAG" = Lys
  | t == "GAU" = Asp
  | t == "GAC" = Asp
  | t == "GAA" = Glu
  | t == "GAG" = Glu
  | t == "UGU" = Cys
  | t == "UGC" = Cys
  | t == "UGA" = Opal
  | t == "UGG" = Trp
  | t == "CGU" = Arg
  | t == "CGC" = Arg
  | t == "CGA" = Arg
  | t == "CGG" = Arg
  | t == "AGU" = Ser
  | t == "AGC" = Ser
  | t == "AGA" = Arg
  | t == "AGG" = Arg
  | t == "GGU" = Gly
  | t == "GGC" = Gly
  | t == "GGA" = Gly
  | t == "GGG" = Gly
  | otherwise = error "Invalid Codon"

codonToChar :: RNACodon -> Char
codonToChar r
 | r == Ala   = 'A' 
 | r == Arg   = 'R' 
 | r == Asn   = 'N' 
 | r == Asp   = 'D' 
 | r == Gly   = 'G' 
 | r == Ser   = 'S' 
 | r == Trp   = 'W' 
 | r == Cys   = 'C' 
 | r == Tyr   = 'Y' 
 | r == His   = 'H' 
 | r == Gln   = 'Q' 
 | r == Lys   = 'K' 
 | r == Glu   = 'E' 
 | r == Pro   = 'P' 
 | r == Thr   = 'T' 
 | r == Val   = 'V' 
 | r == Met   = 'M' 
 | r == Ile   = 'I' 
 | r == Leu   = 'L' 
 | r == Phe   = 'F' 
 | r == Opal  = ' ' 
 | r == Ochre = ' ' 
 | r == Amber = ' ' 
 | otherwise = error "Invalid Codon"

charToCodon :: Char -> RNACodon
charToCodon r
 | r == 'A' = Ala    
 | r == 'R' = Arg    
 | r == 'N' = Asn    
 | r == 'D' = Asp    
 | r == 'G' = Gly    
 | r == 'S' = Ser    
 | r == 'W' = Trp    
 | r == 'C' = Cys    
 | r == 'Y' = Tyr    
 | r == 'H' = His    
 | r == 'Q' = Gln    
 | r == 'K' = Lys    
 | r == 'E' = Glu    
 | r == 'P' = Pro    
 | r == 'T' = Thr    
 | r == 'V' = Val    
 | r == 'M' = Met    
 | r == 'I' = Ile    
 | r == 'L' = Leu    
 | r == 'F' = Phe    
 | r == ' ' = Opal   
 | r == ' ' = Ochre  
 | r == ' ' = Amber  
 | otherwise = error "Invalid Codon"


