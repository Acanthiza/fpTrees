  
# Create a colour palette for n groups

  col_pal <-  function(n) {
    if (n <= 8) {
      RColorBrewer::brewer.pal(n, "Set2")
    } else {
      hcl(h=seq(0,(n-1)/(n),length=n)*360,c=100,l=65,fixup=TRUE)
    }
  }
  

# turn a vector into a comma separated list of values with a penultimate 'and'
  vec_to_sentence <- function(x,sep=",") {
    
    x[!is.na(x)] %>%
      paste(collapse = "JOINSRUS") %>%
      (function(x) if(sep == ";") {
        
          stringi::stri_replace_last_regex(x,"JOINSRUS", paste0(sep," and ")) %>%
            str_replace_all("JOINSRUS",paste0(sep," "))
      
      } else {
        
        stringi::stri_replace_last_regex(x,"JOINSRUS", " and ") %>%
          str_replace_all("JOINSRUS",paste0(sep," "))
        
        }
      )

  }
  

# https://github.com/ateucher/useful_code/blob/master/R/numbers2words.r

  numbers2words <- function(x){
    ## Function by John Fox found here:
    ## http://tolstoy.newcastle.edu.au/R/help/05/04/2715.html
    ## Tweaks by AJH to add commas and "and"
    helper <- function(x){

      digits <- rev(strsplit(as.character(x), "")[[1]])
      nDigits <- length(digits)
      if (nDigits == 1) as.vector(ones[digits])
      else if (nDigits == 2)
        if (x <= 19) as.vector(teens[digits[1]])
      else trim(paste(tens[digits[2]],
                      Recall(as.numeric(digits[1]))))
      else if (nDigits == 3) trim(paste(ones[digits[3]], "hundred and",
                                        Recall(makeNumber(digits[2:1]))))
      else {
        nSuffix <- ((nDigits + 2) %/% 3) - 1
        if (nSuffix > length(suffixes)) stop(paste(x, "is too large!"))
        trim(paste(Recall(makeNumber(digits[
          nDigits:(3*nSuffix + 1)])),
          suffixes[nSuffix],"," ,
          Recall(makeNumber(digits[(3*nSuffix):1]))))
      }
    }
    trim <- function(text){
      #Tidy leading/trailing whitespace, space before comma
      text=gsub("^\ ", "", gsub("\ *$", "", gsub("\ ,",",",text)))
      #Clear any trailing " and"
      text=gsub(" and$","",text)
      #Clear any trailing comma
      gsub("\ *,$","",text)
    }
    makeNumber <- function(...) as.numeric(paste(..., collapse=""))
    #Disable scientific notation
    opts <- options(scipen=100)
    on.exit(options(opts))
    ones <- c("", "one", "two", "three", "four", "five", "six", "seven",
              "eight", "nine")
    names(ones) <- 0:9
    teens <- c("ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen",
               "sixteen", " seventeen", "eighteen", "nineteen")
    names(teens) <- 0:9
    tens <- c("twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty",
              "ninety")
    names(tens) <- 2:9
    x <- round(x)
    suffixes <- c("thousand", "million", "billion", "trillion")
    if (length(x) > 1) return(trim(sapply(x, helper)))
    helper(x)
  }

# Fix html widget when not displayed
  widgetFix <- function(inputFile,outputFile){
    a = readLines(inputFile)
    output = paste(a,collapse="\n")
    output = gsub(">\n\n</div>","></div>",output)
    writeLines(output,outputFile)
    invisible(NULL)
  }
