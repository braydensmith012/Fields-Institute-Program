import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class GeneralSDEMathematica {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		//EDITED FOR SYMMETRY, SYMMETRIC UNDER FLIPS TOGGLE
		
		// MAIN PARAMETERS HERE FOR MODEL
		
		
		// The number of matrices involved in the model (A, B, C, etc.)
		// E.g. for the ABAB model this would be 2
		int lNum = 2;
		
		// Extra words of the model go here (e.g. ABAB, AABB)
		// Quartic terms are automatically added
		// E.g. for the augmented ABAB model, you would add "ABAB" and "AABB"
		String[] extraWords =  new String[] {"ABAB", "AABB"};
		
		// this variable is described by it's name. It generates all SDE's made with
		// words of length less than or equal to this number.
		// 5 or less: pretty much instant in mathematica
		// 7 or less: under a minute in mathematica
		// 9 or less: can possibly take very long
		int maxLengthofWordsGeneratingSDEs = (2 * 2)+1;
		
		// Set to true if the model is symmetric under any switch of letters/matrices
		// in the symmetry group on the letters/matrices
		// If false, the model is assumed to be symmetric only on cyclic switches of letters
		// E.g. A -> B -> C -> A
		// Note that the two options are equivalent if only two or less matrices are involved
		// in the model
		boolean fullSymmetry = true;
		
		
		// Flip symmetry
		boolean flipSymmetry = true;
		
		
		// filename of output file
		String filename = "output.txt";
		
		// END OF MAIN PARAMETERS
		
		
		
		
		
		
		
		
		
		String[][][] firstSDE;
        String[][][] secondSDE;
        String equations = "";
        String Word;
        String printLine = "";
        String printLineX = "";
        int number;
        String[] allWords;
        
        
        // Gets highest X variables needed
        int maximum = 3;
        
        for (int i = 0; i < extraWords.length; i++) {
        	if (extraWords[i].length() > maximum) {
        		maximum = extraWords[i].length();
        	}
        }
        
        
        
        
        maximum += maxLengthofWordsGeneratingSDEs;
        
        
        int[][] extraW = new int[extraWords.length][];
		
		for (int i = 0; i < extraWords.length; i++) {
			extraW[i] = convertLettersToNumbers(extraWords[i], lNum);
		}
				
		
        
        
        
		
		LinkedHashSet<String> set = new LinkedHashSet<String>();
		String[] evenABStrings = null;
        for (int length = 1; length < maximum + 1; length++) {
        	evenABStrings = generateEvenABStrings(length, lNum).toArray(new String[0]);
        
        	System.out.println(length);
        
       for (int i = 0; i < evenABStrings.length; i++) {
    	   if (fullSymmetry) {
    		   evenABStrings[i] = firstWordUltraSymmetry(evenABStrings[i], convertLettersToNumbers(evenABStrings[i], lNum), lNum, flipSymmetry);
    	   } else {
    		   evenABStrings[i] = firstWordSymmetry(evenABStrings[i], convertLettersToNumbers(evenABStrings[i], lNum), lNum, flipSymmetry);
    	   }
    	    set.add(evenABStrings[i]);
       	}
       
        }
        // End of getting highest X variables needed
        
        
        
        
		try (PrintWriter writer = new PrintWriter(new FileWriter(filename))) {
				int[] integers;	
				
				
				
				String[] varX = new String[set.size()];
				for (int i = 1; i < set.size() + 1; i++) {
					varX[i - 1] = "x" + i;
					System.out.println("x" + i + " = " + set.toArray()[i-1]);
				}
				
				writer.print("{");
				for (int i = 0; i < set.size(); i++) {
					writer.print(varX[i] + ", ");
				}
				
				
				writer.println("g}");
				
				

				
				for (int i = 0; i < extraWords.length; i++) {
					writer.println("t" + (i + 1) + " = g;");
				}
				writer.println("t" + (extraWords.length + 1) + " = g;");
				
				
				
				
				
				writer.print("equations = {");
				
				
				
				
				
				for (int max = 0; max < maxLengthofWordsGeneratingSDEs + 1; max++) {
				
				
				
				number = 0;
				
				
				
				
				allWords = generateWords(lNum, max);
				
				
				
				integers = new int[max + 1];
		        
		        
				
				
				// goes over each word
		        for (int i = 0; i < allWords.length; i++) {
		        	
		        	
		        	
		        	
		        	if (checkOdd(allWords[i], lNum)) {
		        		Word = allWords[i];
		        		
		        		// Gets word in number form as well
		        		integers = convertLettersToNumbers(Word, lNum);
		        		

		        		// Word = convertNumbersToLetters(integers);
		        	
		        	    
		        		
		        		// Gets SDE's, removing empty terms	
		        		firstSDE = firstSDE(Word, integers, lNum);
		        		firstSDE = firstToZero(firstSDE, lNum, fullSymmetry, flipSymmetry);
		        
		       			secondSDE = secondSDE(Word, integers, extraWords, extraW, lNum);
		        		secondSDE = secondToZero(secondSDE, lNum, fullSymmetry, flipSymmetry);
		        
		        		// Removes completely empty SDE parts

		        		System.out.println("Word " + Word);
		        		
		        			printLine = printLine + printFirstArray(firstSDE);
		        			printLineX = printLineX + printFirstArrayX(firstSDE, set);
		        		if (areAll3DEmptyStrings(secondSDE) == false) {
		        			printLine = printLine + printSecondArray(secondSDE);
		        			printLineX = printLineX + printSecondArrayX(secondSDE, set);
		        		}
		        		
		        		// Converts SDE's to other code, writes it to file
		        		if (printLine != "") {
		        			if (Word.equals("A")) {
		        				writer.print(printLineX + " == 0");

			        			System.out.println(printLine);
			        			
		        				

				        		equations += "eq" + number + "a" + max + ", ";
				        		
				        		
				        		number++;
		        			} else {
		        			System.out.println(printLine);
		        			

			        		//printing equation
			        		writer.print(", " + printLineX + " == 0");

			        		equations += "eq" + number + "a" + max + ", ";
			        		
			        		
			        		number++;
		        			}
		        		}
		        		
		        		
		        		
		        		
		        		
		        		// Resets variables for next word
		        		printLine = "";
		        		printLineX = "";
		        		integers = new int[max + 1];
		        		
		        	}
		        	
		        }

				} // end of for statement
				
				
				
				
				writer.println("}");
				
				writer.print("solutions = Solve[equations, ");
			    writer.print("{");
			    for (int i = 0; i < varX.length - 1; i++) {
			    	writer.print(varX[i] + ", ");
			    }
			    writer.print(varX[varX.length - 1]);
				writer.println("}]");

				
			    writer.println("rules = solutions[[1]];");
			    writer.println("lhsVars = rules[[All, 1]];");
			    writer.println("rhsVars = Union@Cases[rules[[All, 2]], _Symbol, Infinity];");
			    writer.println("independentVars = Complement[rhsVars, lhsVars]");
				
				
				
				
				
				
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
		        
		
		
		
		// Idea: get in values, set V(H), calculate form of SDE's, use form to solve for all SDE's
		
		// I will need a few things:
		// 1. A set of tracial formula calculators XOR a general formula
		// 2. A way of dynamically adjusting the moments wrt different letter amounts (lexicog?)
		// 3. A way of taking in higher coefficients (t_1, t_2, ..., t_n)
		// 4. Does -N multiply everything in previous model? Need this info. Check binder.
		// 5. Need way of not repeating calculations by converting words into unique standard form
		// 6. How are equations formatted in Mathematica?
		
		// general form of SDE should look like this:
		// 0 = traces multiplied in pairs from inputted word terms +
		// - Tr(WA) + g Tr(WA^3) + 
		// t1 (Tr(W W1(-1 from first A coeff) + W W1 (-1 from second A coeff) + etc) 
		// + ... + tn Tr(W Wn(-1 from A term))
		
		
		// ONE
		// First term of formula should be related to differential of word (function for this)
		// Second term must be in terms of t1,t2, ... tn
		// Second term must be related to differential of V(letters)
		// Perhaps hardcode in quartic terms effect on the SDE's? No need for equations to look nice
		// Use dynamic array to store terms
		
		// TWO
		// Lexicographically sort the moments (x1, x2, etc.)
		// Look at old code for good variable names
		// Have three forms:
		// 1. Number form in an integer array (3,1,2)
		// 2. Letter form in a string (AAABAA)
		// 3. X form in it's place (x5)
		// We work with 1 and 2 in the code, but when printing proper convert to 3
		
		// THREE
		// By taking in text get V(A,B,C, etc.) as variable to be used in ONE
		// Ask for number of extra terms (past quartic)
		// Ask for words of each extra term
		// Ask for max length of word to find SDE for
		// Higher coefficients are only relevant when printing proper at the end
		
		// FOUR
		// Yes, it multiplies everything. Also, the square term is +, every term after -

		// FIVE
		// Check old code on how this was done (try all equivalent forms, get lexicog. first one)
		// Generalize this to multiple letters
		// General symmetry under 
		// 1. letter switch rotations (A to B to C to A)
		// 2. word rotations (ABCD to BCDA)
		// This should work for AABBCC model
		
		// SIX
		// Look at Andrei's python code output for SDE's
		// The required form should be transparent
		
	}
	

    //////////////////////////////////////////////////////////////////
    	
	// Gets list of even moments
    public static List<String> generateEvenABStrings(int length, int lNum) {
        List<String> result = new ArrayList<>();
        int[] count = new int[lNum];
        
        for (int i = 0; i < count.length; i++) {
        	count[i] = 0;
        }
        
        generateEvenABStringsHelper("", count, length, result);
        return result;
    }

    private static void generateEvenABStringsHelper(String currentString, int[] count, int length, List<String> result) {
        String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    	// Base case: if currentString reaches desired length
        if (currentString.length() == length) {
            // Check if counts of A's and B's are even
        	int counter = 0;
        	for (int i = 0; i < count.length; i++) {
        		if (count[i] % 2 == 1) {
        			counter = 1;
        		}
        	}
        	
            if (counter == 0) {
                result.add(currentString);
            }
            return;
        }

        // Recursively generate strings with even A's and B's
        for (int i = 0; i < count.length; i++) {
        	count[i]++;
        	generateEvenABStringsHelper(currentString + alphabet.charAt(i), count, length, result);
        	count[i]--;
        }
        
    }
	

    //////////////////////////////////////////////////////////////////
    
	// Generates strings / words for first part of SDE
	public static String[][][] firstSDE(String Word, int[] w, int lNum)  {
			
			String[][][] words = new String[(w.length) / lNum][][];
			
			// Creates storage based on A coefficients by only looking at even indices
			for (int y=0; y < (w.length)/lNum; y++) {
				words[y] = new String[w[(lNum*(y))]][2];
			}
			
			
			int r = 0;
			
			// Only goes over A coefficients by only looking at even indices
			for (int i = 0; i < (w.length)/lNum; i++) {
				for (int k = 0; k < (lNum*i); k++) {
					r += w[k];
				}
				for (int j = 0; j < w[(lNum*i)]; j++) {
					words[i][j][0] = Word.substring(0, r + j);
					words[i][j][1] = Word.substring(r + j+1);
				}
				
				r = 0;
			}
			
			// For words, the first index represents the term group, the second the term, and
			// the third the tracial part of the term
			return words;
		}

	// Generates strings / words for second part of SDE
	public static String[][][] secondSDE(String Word, int[] w, String[] extraWords, int[][] extraW, int lNum) {
		
		// First index represents the extra word, the second the section of the word
		// and the third the terms of the word
		String[][][] words = new String[2 + extraWords.length][][];
		
		words[0] = new String[1][1];
		words[1] = new String[1][1];
		
		words[0][0][0] = Word + "A";
		words[1][0][0] = Word + "AAA";
 
		
		String wordHolder;
		String wordHolder0;
		int holder;
		
		// 1. For every extraWord, use the number form of the extraWord to get the a coefficient
		// number counter to define the second index
		// 2. Use the orders of the A terms of extraWord to define the final index
		for (int i = 0; i < extraWords.length; i++) {
			words[i + 2] = new String[extraW[i].length / lNum][];
			for (int j = 0; j < extraW[i].length / lNum; j++) { 
				words[i+2][j] = new String[extraW[i][lNum*j]];
				// 2 i,j orders print for both
			}
		}
		
		

		
		// 1. first, check if this section in word is empty
		// 2. If not empty, then subtract one and rotate the rest, appending it to Word
		// 3. If it is empty, skip the section
		// 4. Go to the next section
		// 5. Go to the next word
		int counter = 0;
		for (int i = 0; i < extraWords.length; i++) {
			for (int j = 0; j < extraW[i].length / lNum; j++) { 
				holder = 0;
				counter = 0;
				
				// finds all orders upto the term and sums them into holder
				for (int m = 0; m < j * lNum; m++) {
					holder += extraW[i][m];
				}
				for (int k = 0; k < extraW[i].length; k++) {
					counter += extraW[i][k];
				}
				
				// shifts word over by that number and stores it
				wordHolder = cyclicShift(extraWords[i], counter - holder);
				
				// gets everything past first letter
				wordHolder = wordHolder.substring(1);
				
				// shifted word without first A
				wordHolder0 = wordHolder;
				
				// repeats the number of the order of the jth A at word i
				for (int k = 0; k < words[i+2][j].length; k++) {
					
					//stores everything past first k-1 letters of shifted word without an A
					wordHolder = wordHolder.substring(k);
					
					// Adds k A's to end
					for (int n = 0; n < k; n++) {
						wordHolder = wordHolder + "A";
					}
					
					words[i+2][j][k] = Word + wordHolder;
					wordHolder = wordHolder0;
				}
			}
		}
		
		return words;
	}
	
    //////////////////////////////////////////////////////////////////
    
    // Cyclically shifts Word by k
    public static String cyclicShift(String str, int k) {
        if (str == null || str.length() == 0 || k == 0) {
            return str;
        }

        int len = str.length();
        k = k % len; // Ensure the shift is within the string length

        if (k < 0) {
            k += len; // Handle negative shift by converting it to positive
        }

        // Perform the cyclic shift
        String shiftedStr = str.substring(len - k) + str.substring(0, len - k);
        return shiftedStr;
    }
    
    // Converts Word into it's lexicographically smallest form
    public static String firstWordSymmetry(String Word, int[] w, int lNum, boolean flipSymmetry) {

    	String str1 = "";
    	String str2 = "";
    	String str3 = "";
    	StringBuilder reversed;
    	int[] orders = new int[lNum];
    	int[] wNew = new int[w.length + lNum];
    	
    	
    	// Need to update for possibility of different types of symmetry
    	// e.g. full symmetry group symmetry
    	// something like:
    	List<String> symmetry;
    	
    	symmetry = permuteFirstKLetters(Word, lNum);
    	
    	
    	
    	
    	// now we check if the transpose is smaller or not
    	// idea: take word.length, and use it to get array of all
    	// transposes and cyclic shifts, then sort and get smallest 
    	String holder = Word;
    	for (String s : symmetry) {
    		// if parity odd, flips
    		/*
    		if (Integer.valueOf(s.substring(s.indexOf('!')+1)) % 2 == 1) {
    			s = s.substring(0, s.indexOf('!'));
    			reversed = new StringBuilder(s);
        		s = reversed.reverse().toString();
    		} else {
    			s = s.substring(0, s.indexOf('!'));
    		}
    		*/
    		s = s.substring(0, s.indexOf('!'));
    		
    	for (int i = 0; i < Word.length(); i++) {
    		
    		//cyclic shift goes here, stored into str1
    		str1 = cyclicShift(s, i);
    		//transposition goes here, stored into str2
    		str2 = str1;
    		if (flipSymmetry) {
    		reversed = new StringBuilder(str1);
    		str2 = reversed.reverse().toString();
    		//comparison goes here between two
    		if (isSmaller(str1, str2)) {
        		str2 = str1;
        	}
    		}
    		//comparison with smallest word goes here (str3)
    		if (i == 0 || isSmaller(str2, str3)) {
    			str3 = str2;
    		}
    	}
    	if (isSmaller(str3, holder)) {
    		holder = str3;
    	}
    	}
    	
    	
    	
    	return holder;
    }
    
    // Full symmetric group of switches of letters
    public static String firstWordUltraSymmetry(String Word, int[] w, int lNum, boolean flipSymmetry) {
    	
    	String str1 = "";
    	String str2 = "";
    	String str3 = "";
    	StringBuilder reversed;
    	int[] orders = new int[lNum];
    	int[] wNew = new int[w.length + lNum];
    	
    	
    	// Need to update for possibility of different types of symmetry
    	// e.g. full symmetry group symmetry
    	// something like:
    	List<String> symmetry;
    	
    	symmetry = permuteFirstKLetters(Word, lNum);
    	
    	
    	
    	
    	// now we check if the transpose is smaller or not
    	// idea: take word.length, and use it to get array of all
    	// transposes and cyclic shifts, then sort and get smallest 
    	String holder = Word;
    	for (String s : symmetry) {
    		s = s.substring(0, s.indexOf('!'));
    	for (int i = 0; i < Word.length(); i++) {
    		
    		//cyclic shift goes here, stored into str1
    		str1 = cyclicShift(s, i);
    		//transposition goes here, stored into str2
    		str2 = str1;
    		if (flipSymmetry) {
    		reversed = new StringBuilder(str1);
    		str2 = reversed.reverse().toString();
    		//comparison goes here between two
    		if (isSmaller(str1, str2)) {
        		str2 = str1;
        	}
    		}
    		//comparison with smallest word goes here (str3)
    		if (i == 0 || isSmaller(str2, str3)) {
    			str3 = str2;
    		}
    	}
    	if (isSmaller(str3, holder)) {
    		holder = str3;
    	}
    	}
    	
    	
    	
    	return holder;
    }

    // Finds which of two strings comes lexicographically first
    public static boolean isSmaller(String str1, String str2) {
    	
    	int len = Math.min(str1.length(), str2.length());
    	if (str1.length() < str2.length()) {
    		return true;
   		}
    	for (int j = 0; j < len; j++) {
    	    if (str1.charAt(j) < str2.charAt(j)) {
    	    	return true;
    	    } else if (str1.charAt(j) > str2.charAt(j)) {
    	        return false;
    	    }
    	}
    	return false;     
    }

    ////////////////////////////////////////////////////////////////// 
    
    // Function to generate permutations of the first k letters of the alphabet and apply them to the string
    public static List<String> permuteFirstKLetters(String str, int k) {
        // Validate k
    	String perm1;
        if (k > 26 || k < 1) {
            throw new IllegalArgumentException("k must be between 1 and 26");
        }

        // Generate the first k letters of the alphabet
        String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        String firstKLetters = alphabet.substring(0, k);

        // Generate all permutations of these k letters
        List<String> permutations = new ArrayList<>();
        generatePermutations(firstKLetters.toCharArray(), 0, permutations, 0);

        // Apply each permutation to the string
        Set<String> results = new HashSet<>();
        for (String perm : permutations) {
        	
        	perm1 = perm.substring(0, perm.indexOf('!'));
            results.add(applyPermutation(str, perm1, k) + "!" + perm.substring(perm.indexOf('!')+1));
        }

        return new ArrayList<>(results);
    }

    // Recursive function to generate permutations
    private static void generatePermutations(char[] str, int l, List<String> permutations, int count) {
        if (l == str.length - 1) {
            permutations.add(new String(str) + "!" + count);
        } else {
            for (int i = l; i < str.length; i++) {
                swap(str, l, i);
                count++;
                generatePermutations(str, l + 1, permutations, count);
                swap(str, l, i); // Backtrack
            }
        }
    }

    // Helper function to swap characters in a char array
    private static void swap(char[] str, int i, int j) {
        char temp = str[i];
        str[i] = str[j];
        str[j] = temp;
    }

    // Function to apply a permutation to the first k letters of the string
    private static String applyPermutation(String str, String perm, int k) {
        Map<Character, Character> mapping = new HashMap<>();
        String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        for (int i = 0; i < k; i++) {
            mapping.put(alphabet.charAt(i), perm.charAt(i));
        }

        char[] result = str.toCharArray();
        for (int i = 0; i < result.length; i++) {
            if (mapping.containsKey(result[i])) {
                result[i] = mapping.get(result[i]);
            }
        }

        return new String(result);
    }
    
    //////////////////////////////////////////////////////////////////
        
	// Converts letter form to x form (AAABBA -> A3B2A)
    public static String compressString(String input) {
    	if (input == null || input.isEmpty()) {
            return "";
        }

        
        StringBuilder compressed = new StringBuilder();
        char currentChar = input.charAt(0);
        int count = 1;

        for (int i = 1; i < input.length(); i++) {
            char nextChar = input.charAt(i);
            if (nextChar == currentChar) {
                count++;
            } else {
            	if (count > 1) {
            		compressed.append(currentChar).append(count);
            	} else {
            		compressed.append(currentChar);
            	}
                currentChar = nextChar;
                count = 1;
            }
        }

        // Append the last group of characters
        if (count > 1) {
    		compressed.append(currentChar).append(count);
    	} else {
    		compressed.append(currentChar);
    	}
        return compressed.toString();
    }

	// Method to convert numbers into letters (3,1,2 -> AAABCC) 
    public static String convertNumbersToLetters(int[] numbers, int lNum) {
    	 StringBuilder result = new StringBuilder();
    	 String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    	    int currentLetter = 0;
    	    
    	    for (int number : numbers) {
    	        for (int i = 0; i < number; i++) {
    	            result.append(alphabet.charAt(currentLetter));
    	        }
    	     // Toggle between letters
                currentLetter++;
                if (currentLetter == lNum) {
                	currentLetter = 0;
                }
    	    }
    	    
    	    return result.toString();
    }
	
    // Method to convert letters into numbers (AAABAA -> 3,1,2) (note B^3A -> 3,1)
    public static int[] convertLettersToNumbers(String s, int lNum) {
        List<Integer> result = new ArrayList<>();
        String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        
        if (s == null || s.isEmpty()) {
            return new int[0];
        }
        
        char currentChar = s.charAt(0);
        int count = 1;
        int x;
        
        for (int i = 0; i < alphabet.indexOf(currentChar); i++) { 
        	result.add(0);
        }

        for (int i = 1; i < s.length(); i++) {
            char nextChar = s.charAt(i);
            if (nextChar == currentChar) {
                count++;
            } else {
                result.add(count);
                
                x = alphabet.indexOf(currentChar) + 1;
                while (x != -1) {
                	
                    if (x == lNum) {
                    	
                    	x = 0;
                    }
                    
                    if (alphabet.indexOf(nextChar) != x) {
                    	x++;
                    	result.add(0);
                    } else {
                    	x = -1;
                    }
                }
                
                
                
                currentChar = nextChar;
                count = 1;
            }
        }
        
        // Add the last count
        result.add(count);
        
        for (int i = 0; i < (result.size() % lNum); i++) {
        	result.add(0);
        }


        // Convert List<Integer> to int[]
        int[] counts = new int[result.size()];
        for (int i = 0; i < result.size(); i++) {
            counts[i] = result.get(i);
        }

        
        return counts;
    }
    
    // Gets integer place of word
    public static int getMoment(String Word, LinkedHashSet<String> set) {
    	int count = 0;
    	
    	for (String s : set) {
    		if (Word.equals(s)) {
    			return count + 1;
    		}
    		count++;
    	}
    	return -1;
    }
    
    //////////////////////////////////////////////////////////////////    
    
	// checks if word has odd number of A's, B's, C's, etc. or is odd itself
	public static boolean checkOdd(String Word, int lNum) {
	    	int count = 0;
	    	String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	    	
	    	for (int j = 0; j < lNum; j++) {
	    	for (int i = 0; i < Word.length(); i++) {
	            // Check if the current character matches the target character
	            if (Word.charAt(i) == alphabet.charAt(j)) {
	                count++;
	            }
	        }
	    	if (count % 2 == 1) {
	    		return true;
	    	}
	    	count = 0;
	    	}
	    	if (Word.length() % 2 == 0) {
	    		return false;
	    	}
	    	return true;
	    }
		
    //////////////////////////////////////////////////////////////////    
    
    // Iterate through the 3D array and filter out entries containing '0'
    public static String[][][] filter3DArray(String[][][] array) {
        ArrayList<String[][]> filteredList = new ArrayList<>();

        // Iterate through the 3D array and filter out entries containing '0'
        for (String[][] layer : array) {
            ArrayList<String[]> filteredLayer = new ArrayList<>();
            for (String[] row : layer) {
                boolean containsZero = false;
                for (String str : row) {
                    if (str.contains("0")) {
                        containsZero = true;
                        break;
                    }
                }
                if (!containsZero) {
                    filteredLayer.add(row);
                }
            }
            filteredList.add(filteredLayer.toArray(new String[0][]));
        }

        // Convert ArrayList back to 3D array
        String[][][] filteredArray = new String[filteredList.size()][][];
        for (int i = 0; i < filteredList.size(); i++) {
            filteredArray[i] = filteredList.get(i);
        }

        return filteredArray;
    }

    // Iterate through the 3D array and check if any element is not ""
    public static boolean areAll3DEmptyStrings(String[][][] array) {
        
        for (String[][] layer : array) {
            for (String[] row : layer) {
                for (String str : row) {
                    if (!str.equals("")) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    //////////////////////////////////////////////////////////////////
    
    // sees if any 0 terms in first part of SDE, and puts things in lowest terms
    public static String[][][] firstToZero(String[][][] Word, int lNum, boolean fullSymmetry, boolean flipSymmetry) {
    	
    	for (int i = 0; i < Word.length; i++) {
            for (int j = 0; j < Word[i].length; j++) {
               if (checkOdd(Word[i][j][0], lNum) || checkOdd(Word[i][j][1], lNum)) {
            	   Word[i][j][0] = "0";
            	   Word[i][j][1] = "0";
               } else {
            	   if (Word[i][j][0] != "") {
            		   if (fullSymmetry) {
            			   Word[i][j][0] = firstWordUltraSymmetry(Word[i][j][0], convertLettersToNumbers(Word[i][j][0], lNum), lNum,flipSymmetry); 
                		   
            		   } else {
            		   Word[i][j][0] = firstWordSymmetry(Word[i][j][0], convertLettersToNumbers(Word[i][j][0], lNum), lNum, flipSymmetry); 
            		   }}
            	   if (Word[i][j][1] != "") {
            		   if (fullSymmetry) {
            			   Word[i][j][1] = firstWordUltraSymmetry(Word[i][j][1], convertLettersToNumbers(Word[i][j][1], lNum), lNum, flipSymmetry); 
                		   
            		   } else {
            		   Word[i][j][1] = firstWordSymmetry(Word[i][j][1], convertLettersToNumbers(Word[i][j][1], lNum), lNum, flipSymmetry); 
            		   }}
               }
            }
        }
    	
		return filter3DArray(Word);
    }
    
    // sees if any 0 terms in second part of SDE, and puts things in lowest terms
    public static String[][][] secondToZero(String[][][] Word, int lNum, boolean fullSymmetry, boolean flipSymmetry) {
    	for (int i = 0; i < Word.length; i++) {
            for (int j = 0; j < Word[i].length; j++) {
            	for (int k = 0; k < Word[i][j].length; k++) {
            		if (checkOdd(Word[i][j][k], lNum)) {
            			Word[i][j][k] = "0";
                    } else {
                    	if (fullSymmetry) {
                    		Word[i][j][k] = firstWordUltraSymmetry(Word[i][j][k], convertLettersToNumbers(Word[i][j][k], lNum), lNum, flipSymmetry);
                            
                    	} else {
                    		Word[i][j][k] = firstWordSymmetry(Word[i][j][k], convertLettersToNumbers(Word[i][j][k], lNum), lNum, flipSymmetry);
                            
                    	}
                    	}
            	}
            }
    		
        }
    	
		return filter3DArray(Word);
    }

    //////////////////////////////////////////////////////////////////     
    
    // Makes the array of Words
    public static String[] generateWords(int lNum, int n) {
        int totalWords = (int) Math.pow(lNum, n);
        String[] resultArray = new String[totalWords];
        char[] current = new char[n];
        int[] index = {0}; // Use an array to hold the index reference

        generateWordsHelper(lNum, n, 0, current, resultArray, index);
        return resultArray;
    }
    
    // Makes the array of Words
    private static void generateWordsHelper(int lNum, int n, int depth, char[] current, String[] result, int[] index) {
        String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    	if (depth == n) {
            result[index[0]++] = new String(current);
            return;
        }

        for (int i = 0; i < lNum; i++) {
        	current[depth] = alphabet.charAt(i);
        	generateWordsHelper(lNum, n, depth + 1, current, result, index);
        }
    }
    
    //////////////////////////////////////////////////////////////////     
    
    // Function to print a 3D array, putting letter to number-letter form
    public static String printFirstArray(String[][][] array) {
    	String linePrint = "";
    	
        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[i].length; j++) {
                for (int k = 0; k < array[i][j].length; k++) {
                	if (array[i][j][k] != "" && k == 1) {
                		linePrint = linePrint + "(" + compressString(array[i][j][k]) + ")";
                	} else if (array[i][j][k] != "" && array[i][j][1] == "" && k == 0) {
                		linePrint = linePrint + "(" + compressString(array[i][j][k]) + ")";
                	} else if (array[i][j][k] != "" && k == 0) {
                		linePrint = linePrint + "(" + compressString(array[i][j][k]) + ")*";
                	} else if (k == 1 && array[i][j][1] == "" && array[i][j][0] == "" ) {
                		linePrint = linePrint + "1";
                	}
                	
                }
                linePrint = linePrint + " + ";
            }
        }
        return linePrint;
    }

    // Function to print an array, putting letter to number-letter form, with tn coefficients
    public static String printSecondArray(String[][][] array) {
    	String linePrint = "";
    	int[] counter2 = new int[array.length];
    	int counter = 0;
    	for (int i = 0; i < array.length; i++) {
    	for (int m = 0; m < array[i].length; m++) {
			counter2[i] += array[i][m].length;
		}
    	}
    	
    	int holder = 0;
    	for (int i = 0; i < array.length; i++) {
    		if (i == 0 || holder < counter2[i]) {
    			holder = counter2[i];
    		}
    	}
    	
    	
        for (int i = 0; i < array.length; i++) {
        	
        	
        	
            for (int j = 0; j < array[i].length; j++) {
            	for (int k = 0; k < array[i][j].length; k++) {
            		if (i == 0) {
            			linePrint = linePrint + " - " + compressString(array[i][j][k]);
            			
            		} else if (i == 1) { 
            			linePrint = linePrint + " - (" +  "t" + counter + " ) * x" + compressString(array[i][j][k]);
            		} else {

            			
            			linePrint = linePrint + " - (" +  "t" + counter + " / " + holder + " ) * x" + compressString(array[i][j][k]);
            			
            		}
                }
            }
            counter++;

        }
        return linePrint;
    }
    
    
    //////////////////////////////////////////////////////////////////     
        
    
    // Function to print a 3D array, putting letter to x form
    public static String printFirstArrayX(String[][][] array, LinkedHashSet<String> set) {
    	String linePrint = "";
    	
        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[i].length; j++) {
                for (int k = 0; k < array[i][j].length; k++) {
                	if (array[i][j][k] != "" && k == 1) {
                		linePrint = linePrint + "( x" + getMoment(array[i][j][k], set) + " )";
                	} else if (array[i][j][k] != "" && array[i][j][1] == "" && k == 0) {
                		linePrint = linePrint + "( x" + getMoment(array[i][j][k], set) + " )";
                	} else if (array[i][j][k] != "" && k == 0) {
                		linePrint = linePrint + "( x" + getMoment(array[i][j][k], set) + " )*";
                	} else if (k == 1 && array[i][j][1] == "" && array[i][j][0] == "" ) {
                		linePrint = linePrint + "1";
                	}
                	
                }
                linePrint = linePrint + " + ";
            }
        }
        return linePrint;
    }

    // Function to print an array, putting letter to x form, with tn coefficients
    public static String printSecondArrayX(String[][][] array, LinkedHashSet<String> set) {
    	String linePrint = "";

    	int[] counter2 = new int[array.length];
    	int counter = 0;
    	for (int i = 0; i < array.length; i++) {
    	for (int m = 0; m < array[i].length; m++) {
			counter2[i] += array[i][m].length;
		}
    	}
    	
    	int holder = 0;
    	for (int i = 0; i < array.length; i++) {
    		if (i == 0 || holder < counter2[i]) {
    			holder = counter2[i];
    		}
    	}
    	
    	
    	
    	
        for (int i = 0; i < array.length; i++) {
        	
        	
        	
            for (int j = 0; j < array[i].length; j++) {
            	for (int k = 0; k < array[i][j].length; k++) {
            		if (i == 0) {
            			linePrint = linePrint + " - x" + getMoment(array[i][j][k], set);
            		} else if (i == 1) { 
            			linePrint = linePrint + " - (" +  "t" + counter + " ) * x" + getMoment(array[i][j][k], set);
            		} else {
            			linePrint = linePrint + " - (" +  "t" + counter + " / " + holder + " ) * x" + getMoment(array[i][j][k], set);
            		}
                }
            }
            counter++;
        }
        return linePrint;
    }
    
    
    
}