pass_gen = function(n, n_numb, n_lower, n_upper, n_spec,
    str_numb = as.character(0:9),
    str_lower = letters,
    str_upper = LETTERS,
    str_spec = c("!", "@", "#", "$", "%", "^", "&",
        "*", "(", ")","-", "_", "=", "+", "<", ">",
        "[", "]", "{", "}", "`", "~")){

    if (missing(n))
        stop("The length n must be specified.")
    
    filler = NULL
    if (missing(n_numb)){
        filler = c(filler, str_numb)
        n_numb = 0
        }

    if (missing(n_lower)){
        filler = c(filler, str_lower)
        n_lower = 0
        }

    if (missing(n_upper)){
        filler = c(filler, str_upper)
        n_upper = 0
        }

    if (missing(n_spec)){
        filler = c(filler, str_spec)
        n_spec = 0
        }

    if (n < n_numb + n_lower + n_upper + n_spec){
        n = n_numb + n_lower + n_upper + n_spec
#       message("Given n is smaller than sum of individual parts. Setting n to the sum.")
        }

    out = NULL
    out = c(out, sample(str_numb, n_numb))
    out = c(out, sample(str_lower, n_lower))
    out = c(out, sample(str_upper, n_upper))
    out = c(out, sample(str_spec, n_spec))
    if (length(filler) > 0)
        out = c(out, sample(filler, n - n_numb - n_lower - n_upper - n_spec))
    out = paste0(sample(out), collapse = "")
    
    return (out)
    }
