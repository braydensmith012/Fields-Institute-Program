import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

public class symmetricModels {

	

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		// EDITED FOR SYMMETRY
		
		// MAIN PARAMETERS HERE FOR MODEL

		// Can handle anything lower than (4,6,4)
		// Can handle anything lower than (3,6,6)
		// Can handle anything lower than (2,8,8)
		
		int lNum = 2;
		
		int maxLengthWords = 6;
		
		int maxNumTerms = 2;
		
		boolean fullSymmetry = true;
		
		// END OF MAIN PARAMETERS
		
		
		   String[] array;
		    array = new String[] {"ABCDABCD", "BACDBACD", "DBCADBCA"};
		
		System.out.println("" + symmetryModel(array, lNum));
		
		
		try (PrintWriter writer = new PrintWriter(new FileWriter("output2.txt"))) {

		LinkedHashSet<String> set = new LinkedHashSet<String>();
		String[] evenABStrings = null;
		
		// Generates even words
		for (int length = 1; length < maxLengthWords + 1; length++) {
			evenABStrings = generateEvenABStrings(length, lNum).toArray(new String[0]);
			System.out.println(length);

			for (int i = 0; i < evenABStrings.length; i++) {
				if (set.contains(firstWordBasicSymmetry(evenABStrings[i], lNum)) == false) {

					set.add(firstWordBasicSymmetry(evenABStrings[i], lNum));
				}
			}

		}
		
		String[] stringArray = set.toArray(new String[0]);
		
		
		// gets list of combinations
		List<List<String>> comb = getCombinations(stringArray, maxNumTerms);
        String[][] combinations = convertListTo2DArray(comb);
		
        
        
        
        List<String> symmetric = new ArrayList<>();
		
        String holder = "";
        
		for (int i = 0; i < combinations.length; i++) {
			if (symmetryModel(combinations[i], lNum)) {
				
				for (int j = 0; j < combinations[i].length; j++) {
					if (combinations[i][j] != null) {
						holder += "\"" + combinations[i][j] + "\", ";
					}
					
				}
				
				symmetric.add(holder);
				holder = "";
				
			}
		}
		for (String s : symmetric) {
			writer.println(s);
		}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	
	
	public static boolean symmetryModel(String[] model, int lNum) {
		
		boolean detector = false;
		List<String> symmetry;
		String hold;
		
		// Checks if any duplicate terms in model
		for (int i = 0; i < model.length; i++) {
			if (model[i] != null) {
			if (model[i].equals("AA") || model[i].equals("BB") || model[i].equals("AAAA") || model[i].equals("BBBB")) {
				return false;
			}
				
			
			for (int j = 0; j < i; j++) {
				if (model[i].equals(model[j])) {
					return false;
				}
			}
			
			
	    	// list of all permutations of a word
	        symmetry = permuteFirstKLetters(model[i], lNum);
	        
	        for (int j = 0; j < symmetry.size(); j++) {
	        	// some switch of letters of the word in lowest form
	        	hold = firstWordBasicSymmetry(symmetry.get(j), lNum);
	        
	        	for (int k = 0; k < model.length; k++) {
	        		if (hold.equals(model[k])) {
	        			detector = true;
	        		}
	        	}
	        	
	        	if (detector == false) {
	        		System.out.println(symmetry.get(j));
	        		return false;
	        	}
	        	detector = false;
	        	
	        }
	        	
			}	
	    	}
	    	
	    	
			
		
		return true;
	}
	
	
	 public static String[][] convertListTo2DArray(List<List<String>> listOfLists) {
	        if (listOfLists == null || listOfLists.isEmpty()) {
	            return new String[0][0]; // Return an empty 2D array if the input list is empty
	        }

	        // Determine the number of rows and columns
	        int numRows = listOfLists.size();
	        int numCols = listOfLists.stream().mapToInt(List::size).max().orElse(0);

	        // Create the 2D array
	        String[][] array = new String[numRows][numCols];

	        // Populate the 2D array
	        for (int i = 0; i < numRows; i++) {
	            List<String> innerList = listOfLists.get(i);
	            for (int j = 0; j < innerList.size(); j++) {
	                array[i][j] = innerList.get(j);
	            }
	        }

	        return array;
	    }
	
	public static <T extends Comparable<T>> List<List<T>> getCombinations(T[] array, int k) {
        List<List<T>> result = new ArrayList<>();
        Arrays.sort(array); // Sort the array to handle duplicates

        // Generate combinations of all lengths from 1 to k
        for (int i = 1; i <= k; i++) {
            backtrack(result, new ArrayList<>(), array, i, 0);
        }

        return result;
    }

    private static <T extends Comparable<T>> void backtrack(List<List<T>> result, List<T> tempList, T[] array, int k, int start) {
        if (tempList.size() == k) {
            result.add(new ArrayList<>(tempList));
            return;
        }

        for (int i = start; i < array.length; i++) {
            if (i > start && array[i].equals(array[i - 1])) continue; // Skip duplicates
            tempList.add(array[i]);
            backtrack(result, tempList, array, k, i + 1);
            tempList.remove(tempList.size() - 1);
        }
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
    public static String firstWordBasicSymmetry(String Word, int lNum) {
    	String str1;
    	String str2;
    	String str3 = null;
    	StringBuilder reversed;
    	
    	

    	
    	
    	
    	
    	
    	// now we check if the transpose is smaller or not
    	// idea: take word.length, and use it to get array of all
    	// transposes and cyclic shifts, then sort and get smallest 

    	for (int i = 0; i < Word.length(); i++) {
    		
    		//cyclic shift goes here, stored into str1
    		str1 = cyclicShift(Word, i);
    		//transposition goes here, stored into str2
    		/*
    		reversed = new StringBuilder(str1);
    		str2 = reversed.reverse().toString();
    		//comparison goes here between two
    		if (isSmaller(str1, str2)) {
        		str2 = str1;
        	}
    		*/
    		
    		//comparison with smallest word goes here (str3)
    		//to get old change from str1 to str2 below
    		if (i == 0 || isSmaller(str1, str3)) {
    			str3 = str1;
    		}
    		
    	}
    	
    	
    	return str3;
    }
    
    
    
    // Converts Word into it's lexicographically smallest form
    public static String firstWordSymmetry(String Word, int[] w, int lNum) {
    	String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    	String str1;
    	String str2;
    	String str3 = null;
    	StringBuilder reversed;
    	int[] orders = new int[lNum];
    	int[] wNew = new int[w.length + lNum];
    	
    	for (int i = 0; i < lNum; i++) {
    		orders[i] = 0;
    	}
    	
    	
    	
    	
    	// idea:
    	// 1. get list of orders for every letter
    	// 2. shift to maximal order letter
    	// 3. try every possibility and store lexicographically soonest
    	// 4. return lexicographically soonest from the maximal order letters
    	
    	// Gets array of orders of every letter
    	for (int i = 0; i < w.length; i++) {
    		orders[i % lNum] += w[i];
    	}
    	
    	
    	int maxIndex = 0;
    	
    	// finds alphabetical index of largest order letter
    	for (int i = 0; i < lNum; i++) {
    		if (orders[maxIndex] < orders[i]) {
    			maxIndex = i;
    		} 
    	}

    	String holder = "";
    	
    	// finds things of same size
    	int count = 0;
    	for (int i = 0; i < lNum; i++) {
    		if (orders[maxIndex] == orders[i]) {
    			holder += alphabet.charAt(i);
    			count++;
    		} 
    	}

    	
    	String[] newWord = new String[count];
    	
    	
    	// gets the shifted version of the word in number form
    	for (int j = 0; j < holder.length(); j++) {
    	
    	for (int i = 0; i < wNew.length; i++) {
    		if (i < lNum - alphabet.indexOf(holder.charAt(j))) {
    			wNew[i] = 0;
    		} else if (i < (lNum - alphabet.indexOf(holder.charAt(j))) + w.length) {
    			wNew[i] = w[i - (lNum - alphabet.indexOf(holder.charAt(j)))];
    		} else {
    			wNew[i] = 0;
    		}
    	}
    	newWord[j] = convertNumbersToLetters(wNew, lNum);
    	}
    	
    	
    	
    	
    	
    	
    	// now we check if the transpose is smaller or not
    	// idea: take word.length, and use it to get array of all
    	// transposes and cyclic shifts, then sort and get smallest 
    	for (int j = 0; j < newWord.length; j++) { 
    		for (int i = 0; i < newWord[j].length(); i++) {
    		
    		//cyclic shift goes here, stored into str1
    		str1 = cyclicShift(newWord[j], i);
    		//transposition goes here, stored into str2
    		reversed = new StringBuilder(str1);
    		str2 = reversed.reverse().toString();
    		//comparison goes here between two
    		if (isSmaller(str1, str2)) {
        		str2 = str1;
        	}
    		
    		//comparison with smallest word goes here (str3)
    		if (i == 0 || isSmaller(str2, str3)) {
    			str3 = str2;
    		}
    		
    		}
    		newWord[j] = str3;
    	}
    	
    	String holder2 = newWord[0];
    	
    	for (int j = 1; j < newWord.length; j++) { 
    		if (isSmaller(newWord[j], holder2)) {
    			holder2 = newWord[j];
    		}
    	}
    	
    	return holder2;
    }
    
    // Full symmetric group of switches of letters
    public static String firstWordUltraSymmetry(String Word, int[] w, int lNum) {
    	
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
    	for (int i = 0; i < Word.length(); i++) {
    		
    		//cyclic shift goes here, stored into str1
    		str1 = cyclicShift(s, i);
    		//transposition goes here, stored into str2
    		reversed = new StringBuilder(str1);
    		str2 = reversed.reverse().toString();
    		//comparison goes here between two
    		if (isSmaller(str1, str2)) {
        		str2 = str1;
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
        if (k > 26 || k < 1) {
            throw new IllegalArgumentException("k must be between 1 and 26");
        }

        // Generate the first k letters of the alphabet
        String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        String firstKLetters = alphabet.substring(0, k);

        // Generate all permutations of these k letters
        List<String> permutations = new ArrayList<>();
        generatePermutations(firstKLetters.toCharArray(), 0, permutations);

        // Apply each permutation to the string
        Set<String> results = new HashSet<>();
        for (String perm : permutations) {
            results.add(applyPermutation(str, perm, k));
        }

        return new ArrayList<>(results);
    }

    // Recursive function to generate permutations
    private static void generatePermutations(char[] str, int l, List<String> permutations) {
        if (l == str.length - 1) {
            permutations.add(new String(str));
        } else {
            for (int i = l; i < str.length; i++) {
                swap(str, l, i);
                generatePermutations(str, l + 1, permutations);
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
    
    //////////////////////////////////////////////////////////////////     
    
	
	
	
	
}
