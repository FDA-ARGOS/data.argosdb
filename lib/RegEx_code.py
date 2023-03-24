#RegEx Expression Patterns and Associated Code 

# Match everything from start to end of the string 
import re

def use_regex(input_text):
    pattern = re.compile(r"^.*$")
    return pattern.match(input_text)
  
# Combination (Alphanumeric Charaters + Characters)
import re

def use_regex(input_text):
    pattern = re.compile(r"^([A-Za-z0-9]+( [A-Za-z0-9]+)+)$", re.IGNORECASE)
    return pattern.match(input_text)

 # Combination (Decimal number + Character)
import re

def use_regex(input_text):
    pattern = re.compile(r"^([0-9]*\.[0-9]+(\|[0-9]*\.[0-9]+)+)$", re.IGNORECASE)
    return pattern.match(input_text)

# Combination (Number + Character) 
import re

def use_regex(input_text):
    pattern = re.compile(r"^([0-9]+(\|[0-9]+)+)$", re.IGNORECASE)
    return pattern.match(input_text)
  
 # Number 
 import re

def use_regex(input_text):
    pattern = re.compile(r"[0-9]+", re.IGNORECASE)
    return pattern.match(input_text)
  
 # Alphanumeric
import re

def use_regex(input_text):
    pattern = re.compile(r"^[A-Za-z0-9]+$", re.IGNORECASE)
    return pattern.match(input_text)
  
 # Multiple Characters 
import re

def use_regex(input_text):
    pattern = re.compile(r"^[A-Za-z]+$", re.IGNORECASE)
    return pattern.match(input_text)

# Decimal Number 
import re

def use_regex(input_text):
    pattern = re.compile(r"^[0-9]*\.[0-9]+$", re.IGNORECASE)
    return pattern.match(input_text)
  
  
