#!/usr/bin/env python3
"""
Peptidase Cleavage Analysis Tool

This script performs two main functions:
1. Literature Review: Searches for cleavage specificity of a given peptidase using Perplexity API.
2. Cleavage Prediction: Finds cleavage sites in a given protein sequence based on a pattern.
"""

import sys
import os
import argparse
import re
import json
from typing import List, Dict, Optional, Tuple, Union, Any

# Add repository root to path to allow imports from other skills if needed
# But we are keeping this skill relatively self-contained now.
current_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(current_dir)))
# sys.path.append(repo_root)

# --- Perplexity Search Implementation (Self-Contained) ---

def check_dependencies():
    """Check if required packages are installed."""
    try:
        import litellm
        return True
    except ImportError:
        return False

def check_api_key() -> Optional[str]:
    """Check if OpenRouter API key is configured."""
    api_key = os.environ.get("OPENROUTER_API_KEY")
    return api_key

def search_with_perplexity(
    query: str,
    model: str = "openrouter/perplexity/sonar-pro",
    max_tokens: int = 4000,
    temperature: float = 0.2,
    verbose: bool = False
) -> Dict[str, Any]:
    """
    Perform a search using Perplexity models via LiteLLM and OpenRouter.
    Dependencies: litellm
    """
    try:
        from litellm import completion
    except ImportError:
        return {
            "success": False,
            "error": "LiteLLM not installed. Run: uv pip install litellm"
        }

    # Check API key
    api_key = check_api_key()
    if not api_key:
        return {
            "success": False,
            "error": "OpenRouter API key not configured. Set OPENROUTER_API_KEY environment variable."
        }

    if verbose:
        print(f"Model: {model}", file=sys.stderr)
        print(f"Query: {query}", file=sys.stderr)

    # Prepend openrouter/ to model name if not already present
    if not model.startswith("openrouter/"):
        model = f"openrouter/perplexity/{model}"

    try:
        # Perform the search using LiteLLM
        response = completion(
            model=model,
            messages=[{
                "role": "user",
                "content": query
            }],
            max_tokens=max_tokens,
            temperature=temperature
        )

        # Extract the response
        result = {
            "success": True,
            "query": query,
            "model": model,
            "answer": response.choices[0].message.content,
            "usage": {
                "prompt_tokens": response.usage.prompt_tokens,
                "completion_tokens": response.usage.completion_tokens,
                "total_tokens": response.usage.total_tokens
            }
        }

        # Check if citations are available in the response
        if hasattr(response.choices[0].message, 'citations'):
            result["citations"] = response.choices[0].message.citations

        return result

    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "query": query,
            "model": model
        }

# --- End Perplexity Search Implementation ---

class PeptidaseAnalyzer:
    def __init__(self, verbose: bool = False):
        self.verbose = verbose

    def literature_review(self, peptidase_name: str) -> Dict[str, str]:
        """
        Search literature for peptidase specificity and cleavage patterns.
        """
        if not check_dependencies():
            return {
                "success": False,
                "error": "LiteLLM not installed. Cannot perform literature review. Please install 'litellm'."
            }

        query = (
            f"What is the specific cleavage site recognition sequence or pattern for the peptidase '{peptidase_name}'? "
            "Please provide the consensus sequence, P1-P1' notation, and any rules about residues that prevent cleavage "
            "(e.g., 'cleaves after Lys or Arg, unless followed by Pro'). "
            "If possible, provide a Regex pattern representing this specificity."
        )

        if self.verbose:
            print(f"[*] Searching literature for {peptidase_name} specificity...")

        try:
            # Using a reasoning model to get better structural understanding if possible, otherwise default
            response = search_with_perplexity(
                query=query,
                model="sonar-pro", # reliable default
                verbose=self.verbose
            )
            
            if response["success"]:
                return {
                    "success": True,
                    "description": response["answer"],
                    "source": "Perplexity/Literature"
                }
            else:
                return {
                    "success": False,
                    "error": response.get("error", "Unknown error during search")
                }
                
        except Exception as e:
            return {"success": False, "error": str(e)}

    def predict_cleavage_sites(self, sequence: str, pattern: str, pattern_type: str = "regex") -> List[Dict]:
        """
        Find cleavage sites in a sequence using a regex pattern.
        """
        sites = []
        
        try:
            # Compile regex
            # We use lookahead/behind assertions often in protease patterns
            regex = re.compile(pattern)
            
            for match in regex.finditer(sequence):
                # Standard convention: Cleavage is after the matched group, 
                # OR if the user provided a complex regex with groups, we might need logic.
                # simpler assumption: The regex matches the P1 residue (or sequence ending at P1).
                
                cleavage_site = match.end()
                
                # Check boundaries
                if 0 < cleavage_site < len(sequence):
                    # Get context (P4-P4')
                    start_ctx = max(0, cleavage_site - 4)
                    end_ctx = min(len(sequence), cleavage_site + 4)
                    context = sequence[start_ctx:end_ctx]
                    
                    # Add marker '|' to context string for visualization
                    # P4 P3 P2 P1 | P1' P2' P3' P4'
                    # We need to insert '|' relative to the context slice
                    rel_pos = cleavage_site - start_ctx
                    context_visual = context[:rel_pos] + "|" + context[rel_pos:]

                    sites.append({
                        "position": cleavage_site, # 1-based index of residue BEFORE cleavage? No, usually 1-based index of P1.
                                                 # But let's stick to 0-based index for 'position is AFTER residue i'
                        "p1_residue": sequence[cleavage_site-1],
                        "p1_prime_residue": sequence[cleavage_site],
                        "context": context_visual,
                        "match_start": match.start(),
                        "match_end": match.end()
                    })
                    
            return sites
            
        except re.error as e:
            print(f"[!] Invalid Regex Pattern: {e}")
            return []

    def extract_pattern_from_text(self, text: str) -> Optional[str]:
        """
        Attempt to extract a regex pattern from the literature text text.
        Very heuristic-based.
        """
        # Look for regex-like patterns in code blocks or quotes
        # This is a simple heuristic and might need user verification
        
        # Check for explicit mention of regex
        regex_matches = re.findall(r"regex:?\s*[`'\"]([^`'\"]+)[`'\"]", text, re.IGNORECASE)
        if regex_matches:
            return regex_matches[0]
            
        return None

def main():
    parser = argparse.ArgumentParser(description="Analyze Peptidase Specificity and Cleavage Sites")
    parser.add_argument("--peptidase", "-p", help="Name of the peptidase (e.g., 'Trypsin', 'Caspase-3')")
    parser.add_argument("--sequence", "-s", help="Target protein sequence to analyze")
    parser.add_argument("--pattern", "-r", help="Regex pattern for cleavage (e.g., '(?<=[KR])(?!P)')")
    parser.add_argument("--mode", choices=["literature_only", "predict_only", "full"], default="full", 
                        help="Mode of operation: literature search only, prediction only (requires pattern), or full end-to-end")
    parser.add_argument("--output", "-o", help="Output file for results (JSON)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output")
    
    args = parser.parse_args()
    
    analyzer = PeptidaseAnalyzer(verbose=args.verbose)
    results = {}
    
    # 1. Literature Review
    if args.mode in ["literature_only", "full"]:
        if not args.peptidase:
            print("[!] Peptidase name required for literature review mode.")
            return
            
        print(f"[*] Analyzing specificity for: {args.peptidase}")
        lit_results = analyzer.literature_review(args.peptidase)
        
        if lit_results["success"]:
            print("\n=== Literature Review Results ===")
            print(lit_results["description"])
            results["literature_review"] = lit_results
            
            # Try to infer pattern if not provided
            if not args.pattern:
                inferred_pattern = analyzer.extract_pattern_from_text(lit_results["description"])
                if inferred_pattern:
                    print(f"\n[*] Inferred Regex Pattern: {inferred_pattern}")
                    results["inferred_pattern"] = inferred_pattern
                    # Use inferred pattern if we are continuing to prediction
                    if args.mode == "full":
                        args.pattern = inferred_pattern
                        
        else:
            print(f"[!] Literature review failed: {lit_results.get('error')}")
            if args.mode == "full" and not args.pattern:
                print("[!] Cannot proceed to prediction without a pattern.")
                return

    # 2. Cleavage Prediction
    if args.mode in ["predict_only", "full"]:
        if not args.sequence:
            if args.mode == "predict_only":
                print("[!] Sequence required for prediction mode.")
                return
            else:
                print("\n[*] No sequence provided, skipping prediction.")
        
        elif not args.pattern:
             print("[!] No cleavage pattern available. Please provide --pattern or ensure literature review finds one.")
             
        else:
            print(f"\n[*] Predicting cleavage sites in sequence (Length: {len(args.sequence)})")
            print(f"[*] Pattern: {args.pattern}")
            
            sites = analyzer.predict_cleavage_sites(args.sequence, args.pattern)
            results["cleavage_sites"] = sites
            results["cleavage_count"] = len(sites)
            
            print(f"\n=== Prediction Results ({len(sites)} sites found) ===")
            print(f"{'Pos':<6} | {'P1':<3} | {'P1\'':<3} | {'Context':<15}")
            print("-" * 35)
            for site in sites:
                print(f"{site['position']:<6} | {site['p1_residue']:<3} | {site['p1_prime_residue']:<3} | {site['context']}")

    # 3. Output
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\n[*] Results saved to {args.output}")

if __name__ == "__main__":
    main()
