### PC Algorithm for Causal Discovery
### Clark Kaminsky (clarkkam@umich.edu), Mohammed Ombadi (ombadi@umich.edu); 2023


# Setup -- import data
# Set max time lag, significance threshold, max combinations


# Set up time lagged variables


# Loop over all variables
    # Identify parents
    # Loop over the max number of parents, check if there are enough parents to consider
        #n = 0, N(Xt) = {Xt-1 ... Xt-taumax, Yt-1 ... Yt-taumax, ...}
        # While |N(Xt)| <= n
            # Select n dimensional subset of N(Xt)
            # Calculate conditional mutual information and partial correlation
    # Remove links below significance threshold