function histogram_cic_2d( q1, q2, w,
    nbins_1::Int, bins_start_1::Float64, bins_end_1::Float64,
    nbins_2::Int, bins_start_2::Float64, bins_end_2::Float64 )
    """
    Return an 2D histogram of the values in `q1` and `q2` weighted by `w`,
    consisting of `nbins_1` bins in the first dimension and `nbins_2` bins
    in the second dimension.
    Contribution to each bins is determined by the
    Cloud-in-Cell weighting scheme.
    Source: 
    ----------
        https://github.com/openPMD/openPMD-viewer/blob/dev/openpmd_viewer/openpmd_timeseries/utilities.py
    ----------
    
    Parameters:
    ----------
    q1/q2 : float, array
            q1/q2 position of the particles
    w: float, array
            weights (charges) of the particles
    nbins_1/nbins_2 : int
            number of bins (vertices) in the q1/q2 direction
    bins_start_1, bins_end_1, bins_start_2, bins_end_2: float
            start/end value in the q1/q2 direction
    ----------
    
    Returns:
    ----------
    A 2D array of size (nbins_1, nbins_2)
    ----------
    """
    # Define various scalars
    bin_spacing_1 = (bins_end_1-bins_start_1)/(nbins_1-1)
    inv_spacing_1 = 1.0/bin_spacing_1
    bin_spacing_2 = (bins_end_2-bins_start_2)/(nbins_2-1)
    inv_spacing_2 = 1.0/bin_spacing_2
    n_ptcl = length(w)

    # Allocate array for histogrammed data
    hist_data = zeros( (nbins_1, nbins_2) )

    # Go through particle array and bin the data
    for i in 1:n_ptcl

        # Calculate the index of lower bin to which this particle contributes
        q1_cell = (q1[i] - bins_start_1) * inv_spacing_1
        q2_cell = (q2[i] - bins_start_2) * inv_spacing_2
        
        #i1_low_bin = int( math.floor( q1_cell ) )
        
        i1_low_bin = convert(Int, floor(q1_cell)) 
        i2_low_bin = convert(Int, floor(q2_cell)) 

        # Calculate corresponding CIC shape and deposit the weight
        S1_low = 1.0 - (q1_cell - i1_low_bin)
        S2_low = 1.0 - (q2_cell - i2_low_bin)
        if (i1_low_bin >= 0) && (i1_low_bin < nbins_1)
            if (i2_low_bin >= 0) && (i2_low_bin < nbins_2)
                hist_data[ i1_low_bin+1, i2_low_bin+1 ] += w[i]*S1_low*S2_low
            end
            if (i2_low_bin+1 >= 0) && (i2_low_bin+1 < nbins_2)
                hist_data[ i1_low_bin+1, i2_low_bin+2 ] += w[i]*S1_low*(1.0-S2_low)
            end
        end
        if (i1_low_bin+1 >= 0) && (i1_low_bin+1 < nbins_1)
            if (i2_low_bin >= 0) && (i2_low_bin < nbins_2)
                hist_data[ i1_low_bin+2, i2_low_bin+1 ] += w[i]*(1.0-S1_low)*S2_low
            end
            if (i2_low_bin+1 >= 0) && (i2_low_bin+1 < nbins_2)
                hist_data[ i1_low_bin+2, i2_low_bin+2 ] += w[i]*(1.0-S1_low)*(1.0-S2_low)
            end
        end
    end

    return( hist_data )
                                
end
