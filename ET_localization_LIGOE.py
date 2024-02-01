import sys
import math
from numpy import *
from scipy import special
#import glue
#from pylal.xlal import tools
#from pylal import inject
import matplotlib.pyplot
matplotlib.use('Agg')
import pylab

# set the plotting specs
pylab.rcParams.update({
    #"text.verticalalignment": "center",
    "lines.markersize": 6,
    "lines.markeredgewidth": 1.5,
    "lines.linewidth": 1.0,
    "font.size": 20,
    "axes.titlesize": 20,
    "axes.labelsize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 20,
})



n = 30  #Size of array
ntrials = 100000

N = n**2     #No of values

longmin = -11
longmax = 45
deltalong = (longmax - longmin)/float(n)
latmin = 35
latmax = 72
deltalat = (latmax - latmin)/float(n)




latitude = zeros(N)
longitude = zeros(N)
    
for i in range(n):
    for j in range(n):
     
        latitude[j+n*i] = latmin + j*deltalat     #acts as y coord
            
        longitude[i+j*n] = longmin + j*deltalong    #acts as x coord





def calc_location_response(longitude, latitude, arms):
    """
    Calculate the location and response for a detector with longitude, latitude in degrees
    The angle gives the orientation of the arms and is in degrees from North to East
    """
    phi = radians(longitude)
    theta = radians(latitude)
    angle = radians(arms)
    r = 6.4e6
    location = r * xyz(phi, theta)
    r_hat = location / linalg.norm(location)
    # Take North, project onto earth's surface...
    e_n = array([0,0,1])
    e_n = e_n - r_hat * inner(e_n, r_hat)
    # normalize
    e_n = e_n / linalg.norm(e_n)
    # and calculate east
    e_e = cross(e_n, r_hat)
    # Calculate arm vectors
    u_y = e_e * sin(angle) + e_n * cos(angle)
    u_x = e_e * sin(angle + pi/2) + e_n * cos(angle + pi/2)
    response = array(1./2 * (outer(u_x, u_x) - outer(u_y, u_y)), dtype=float32)
    return location, response



def detectors(g, bangalore=True):
    location = {}
    response = {}    
    location["H1"], response["H1"] = calc_location_response(-(119+24/60+27.565681/3600), 46+27/60+18.527841/3600, 36)
    location["L1"], response["L1"] = calc_location_response(-(90+46/60+27.265294/3600), 30+33/60+46.419531/3600, 118)
    location["V1"], response["V1"] = calc_location_response(10+30/60+16.1878/3600, 43+37/60+53.0921/3600, -19.4326 )
    

    
    location["ET"], response["ET"] = calc_location_response(latitude[g], longitude[g], -19.4326 )
    # AIGO location/response:
    #
    # The Australian International Gravitational Observatory (AIGO)
    # Located at Gingin (115 degrees 42 minutes and 30 seconds east, 31
    # degrees 21 minutes and 30 seconds south)
    #
    location["A1"], response["A1"] = calc_location_response(115 + 42./60, -31 + 21./60, 0)
    #
    # KAGRA location:
    # Here is the coordinates
    # 36.25 degree N, 136.718 degree E
    # and 19 degrees from North to West.
    location["K1"], response["K1"] = calc_location_response(136.718, 36.25, -61.7)
    # INDIGO 
    if bangalore:
        # Could you pl run us the Localisation Plot once more with
        # a location close to Bangalore that is seismically quiet
        # 14 deg 14' N
        # 76 deg 26' E?
        location["I1"], response["I1"] = calc_location_response(76 + 26./60, 14 + 14./60, 270)
    else:
        # Here is the coordinates
        # location: 74deg  02' 59" E 19deg 05' 47" N 270.0 deg (W)
        location["I1"], response["I1"] = calc_location_response(74 + 3./60, 19 + 6./60, 270)
    # 
    return( location, response )
    

################################
# plot response functions
################################
def plot_response(ifos, psi, D, X):
    
    #Plot the network response function for a list of ifos.
    
    dra, ddec= pi/100.0, pi/100.0
    [ra,dec] = mgrid[-pi:pi+dra*0.9:dra, -pi/2:pi/2+ddec*0.9:dra]
    f_rss = zeros(shape(ra))
    for i in xrange(shape(ra)[0]):
        for j in xrange(shape(ra)[1]):
           for ifo in ifos:
               #D = 16*365 + 68
               
               gmst = 0.0#18.697374558 + 24.06570982441908 * D
               
               gha = gmst - ra     #greenwich mean sidereal time (rad) - right ascension (radians)
                
               cosgha = np.cos(gha)
               singha = np.sin(gha)
               cosdec = np.cos(dec)
               sindec = np.sin(dec)
               cospsi = np.cos(psi)
               sinpsi = np.sin(psi)
               
               
               X[0] = -cospsi * singha - sinpsi * cosgha * sindec
               X[1] = -cospsi * cosgha + sinpsi * singha * sindec
               X[2] = sinpsi * cosdec
               
               Y[0] = sinpsi * singha - cospsi * cosgha * sindec
               Y[1] = sinpsi * cosgha + cospsi * singha * sindec
               Y[2] = cospsi * cosdec
               
               #f_plus = f_cross = 0.0
                
                
               for i in range(0, 3):
                   DX = D[i][0] * X[0] + D[i][1] * X[1] + D[i][2] * X[2]
                   DY = D[i][0] * Y[0] + D[i][1] * Y[1] + D[i][2] * Y[2]
                   f_plus += X[i] * DX - Y[i] * DY
                   f_cross += X[i] * DY + Y[i] * DX
               
           f_rss[i][j] += f_plus**2 + f_cross**2
    f_rss = sqrt(f_rss)
    pylab.figure()
    pylab.axes(projection="hammer")
    pylab.grid()
    pylab.contourf(ra,dec,f_rss)
    pylab.colorbar()



  
################################
# helper function to make pairs
################################
def pairs(ifos):
    pairs = []
    for det1 in ifos:
        for det2 in ifos:
            if det1 < det2:
                pairs.append([det1,det2])
    return(pairs)


################################
# co-ordinate transformations 
################################
def xyz(ra, dec):
    """
    ra, dec -> x,y,z
    """
    x = cos(dec) * cos(ra)
    y = cos(dec) * sin(ra)
    z = sin(dec)
    loc = asarray([x,y,z])
    return(loc)

def radec(x, y, z):
    """
    x,y,z -> ra, dec
    """
    dec = arcsin(z)
    ra = arctan2(y,x)
    return(ra, dec)  
    

################################
# Calculate localizations 
################################
def localization(ifos, sigma_t):
    """
    Calculate the localization matrix for a set of detectors, 
    This is equation (8) in the advanced detector localization paper.
    We also return the eigenvectors and sigmas (1/eval)
    """
    M = zeros([3,3])
    #nifo = len(ifos)
    sigma_fac = 0
    for d in ifos:
        sigma_fac += sigma_t[d]**-2
    for [d1, d2] in pairs(ifos):
        M += 1.0/sigma_fac * \
           outer(location[d1] - location[d2], location[d1] - location[d2]) / \
           ( (3e8)**2 * sigma_t[d1]**2 * sigma_t[d2]**2 )
    return(M)


def sigmasq_t(ifo, D, ra, dec, psi, cosi, response, bandwidth = 100, D_horizon = 360):
    """
    Calculate the timing accuracy sigma^2_t for a given detector given by ifo and
    sensitivity given by bandwidth and D_horizon 
    and a source with parameters D, ra, dec, psi, cosi
    """
    gmst = 0.0
    
    
    gha = gmst - ra     #greenwich mean sidereal time (rad) - right ascension (radians)
    
    cosgha = cos(gha)
    singha = sin(gha)
    cosdec = cos(dec)
    sindec = sin(dec)
    cospsi = cos(psi)
    sinpsi = sin(psi)
    
    X=zeros(3)
    Y=zeros(3)
    X[0] = -cospsi * singha - sinpsi * cosgha * sindec
    X[1] = -cospsi * cosgha + sinpsi * singha * sindec
    X[2] = sinpsi * cosdec
    
    Y[0] = sinpsi * singha - cospsi * cosgha * sindec
    Y[1] = sinpsi * cosgha + cospsi * singha * sindec
    Y[2] = cospsi * cosdec
    
    fplus = fcross = 0.0
    
    Dmatrix = response[ifo]
    for i in range(0, 3):
        DX = Dmatrix[i][0] * X[0] + Dmatrix[i][1] * X[1] + Dmatrix[i][2] * X[2]
        DY = Dmatrix[i][0] * Y[0] + Dmatrix[i][1] * Y[1] + Dmatrix[i][2] * Y[2]
        fplus += X[i] * DX - Y[i] * DY
        fcross += X[i] * DY + Y[i] * DX
    
    f_plus = fplus
    f_cross = fcross###########################################################################################
    snr = 8. * D_horizon / D * sqrt(f_plus**2 * (1 + cosi*cosi)**2 / 4 + f_cross**2 * cosi * cosi)
    # calculate timing based on equation (1) in advanced localization paper
    sigmasq_t = (2 * pi * snr * bandwidth)**(-2)
    return(sigmasq_t, snr)



def loc_source2(ifos, D, ra, dec, psi, cosi,response, bandwidth_dict, range_dict, snr_found = 5.0):
    """
    Localization of a source by a network of detectors (ifos), assuming the source has
    parameters D, ra, dec, psi, cosi
    Detectors have different bandwidths and ranges contained in the relevant dictionaries, 
    There's a common SNR threshold for "finding" a signal in a detector
    """
    # work out the localization factors
    sigmasq = {}
    snrsq = 0
    loc_factor = 0
    num_found = 0
    for ifo in ifos:
        sigmasq[ifo], snr = sigmasq_t(ifo, D, ra, dec, psi, cosi, response, bandwidth_dict[ifo],       range_dict[ifo] * 2.26)
        loc_factor += 1./sigmasq[ifo]
        snrsq += snr**2
        if (snr > snr_found) and (ifo != "H2"): 
            num_found += 1
    # localization matrix from equation 8 in advanced localization paper
    M = zeros([3,3])
    for [d1, d2] in pairs(dets):
        M += outer(location[d1]- location[d2],  location[d1]-location[d2]) /(3e8)**2 \
            / sigmasq[d1] / sigmasq[d2] / loc_factor
    return(M, snrsq, num_found)

################################
# Operations on the localization matrix M
################################
def evec_sigma(M):
    """
    For a 3x3 matrix M, calculate the eigenvalues and vectors.
    sigma is defined as the reciprocal of the eigenvalue
    """
    ev, evec = linalg.eig(M)
    epsilon = 1e-10
    sigma = 1/ sqrt(ev + epsilon) 
    evec = evec[:,sigma.argsort()]
    sigma.sort()
    return (evec, sigma)

def project(M, source):
    """
    Project localization matrix to zero out components in direction of source
    This is implementing equations 10 and 11 from the advanced localization paper
    """
    P = identity(3) - outer(source, source)
    PMP = inner(inner(P,M),P)
    return(PMP)

def areas(M):
    """
    Calculate the localization areas for a given set of sigma_t
    Check how many are non-zero and proceed accordingly
    Use results from secion 2.2 of advanced localization paper
    """
    evec, sigma = evec_sigma(M)
    if isnan(sigma[0]) or isinf(sigma[0]):
        # then no localization at all:
        a_best = Inf
        a_worst = Inf
    elif isnan(sigma[1]) or isinf(sigma[1]) or \
      sigma[1] > 100 * sigma[0]:
        # only a single baseline (use eqn 15 in advanced localization paper)
        a_best = a_worst = special.erfinv(0.90)*sqrt(2) * 4 * pi * sigma[0] 
        a_worst = Inf
    else:
        # equation 13 in advanced localization paper
        a_best = 2*pi*2.3 * sigma[0] * sigma[1]
        # equation 14 in advanced localization paper -- as bad as it gets with directions 0&1
        a_worst_2 = 2*pi*2.3* sigma[1] *  sqrt(2*sigma[0])
        # localization with directions 1&2
        a_worst_3 = 2*math.pi*2.3* sigma[1] * sigma[2]
        a_worst = min(a_worst_2, a_worst_3)
    # return numbers in degrees
    return(a_best * (180/pi)**2, a_worst*(180/pi)**2)

################################
# plot the 3+ detector ellipse
################################
def plot_projection(M, source):
    # set up co-ordinates
    r = {}
    dangle=  pi/50.0
    angle= arange(0,2*pi+dangle*0.5,dangle)
    # calculate eigendirections
    evec,sigma = evec_sigma(M)
    if sigma[2] < 1:
        # non-degenerate, project onto the sphere
        PMP = project(M, source)
        evec, sigma = evec_sigma(PMP) 
    # calculate the basis
    x_net = evec[:,0]
    y_net = evec[:,1]
    z_net = evec[:,2]
    # normalization to get 90% area:
    x = inner(source,x_net) + \
       sqrt(2 * log(10)) *  sigma[0] * cos(angle) 
    y = inner(source,y_net) + \
      sqrt(2 * log(10)) *  sigma[1] * sin(angle) 
    z = sqrt(1 - x**2 - y**2) * sign(inner(source,z_net))
    # check that we're not going outside of unit circle:
    bad = x**2 + y**2 > 1
    if sum(bad) > 0:
        x = concatenate((x[bad.argmax():], x[0:bad.argmax()]))
        y = concatenate((y[bad.argmax():], y[0:bad.argmax()]))
        z = concatenate((z[bad.argmax():], z[0:bad.argmax()]))
        bad = concatenate((bad[bad.argmax():], bad[0:bad.argmax()]))
        x = x[~bad]
        y = y[~bad]
        z = z[~bad]
        x = append(concatenate((x, x[::-1])), x[0])
        y = append(concatenate((y, y[::-1])), y[0])
        z = append(concatenate((z, -z[::-1])),z[0])
    for i in xrange(3):
        r[i] = x * x_net[i] + y * y_net[i] + z * z_net[i]
    theta = arcsin(r[2] / sqrt(r[0]**2 + r[1]**2 + r[2]**2) )
    phi = arctan2(r[1],r[0])
    pylab.plot(phi, theta, 'b')

###############

def choices(vals, n): 
    
    if n == len(vals): 
            yield tuple(vals) 
    elif n > 1: 
            n -= 1 
            for i, v in enumerate(vals[:-n]): 
                    v = (v,) 
                    for c in choices(vals[i+1:], n): 
                            yield v + c 
    elif n == 1: 
            for v in vals: 
                    yield (v,) 
    elif n == 0: 
            yield () 
    else: 
            # n < 0 
            raise ValueError(n)

#######################################################################
verbose = False
print_net = True 

################################
## Set up detectors and networks
################################

# Notes: 
# 1) Don't have noise curves for Virgo for early and mid runs.  
#    I used the late (no SRM, ~100 Mpc) curve rescaled upwards by factors
#    of 5, 2 to give ranges of ~20, ~50 Mpc and bandwidth of 58.2 Hz 
# 2) For the early, mid, late runs use the "average" noise spectrum 
#    range and bandwidth for localisation.
# 3) Not sure what the numbers for the "steve" row are here -- some 
#    sort of test case.

bandwidth_dict_all = {
    "early" : {'H1' : 123.7, 'L1' : 123.7, 'V1':  58.2, "I1" : 0 }, 
    "mid" :   {'H1' : 103.7, 'L1' : 103.7, 'V1':  58.2, "I1" : 0 }, 
    "late" :  {'H1' :  96.4, 'L1' :  96.4, 'V1':  58.2, "I1" : 0 }, 
    "final" : {'H1' : 117.4, 'L1' : 117.4, 'V1': 148.9, "I1" : 0 }, 
    "india" : {'H1' : 117.4, 'L1' : 117.4, 'V1': 148.9, "I1" : 117.4 }, 
    "steve" : {'H1' : 100.0, 'L1' : 100.0, 'V1': 100.0, "I1" : 100.0 }, 
    "ET" : {'H1' : 57.9535, 'L1' : 57.9535, 'ET': 97.6393 }
}

range_dict_all = {
    "early" : {'H1' :  60.6, 'L1' :  60.6, 'V1':  21.8, "I1" : 0 }, 
    "mid" :   {'H1' : 100.1, 'L1' : 100.1, 'V1':  54.5, "I1" : 0 }, 
    "late" :  {'H1' : 140.8, 'L1' : 140.8, 'V1': 109.0, "I1" : 0 }, 
    "final" : {'H1' : 197.5, 'L1' : 197.5, 'V1': 128.3, "I1" : 0 }, 
    "india" : {'H1' : 197.5, 'L1' : 197.5, 'V1': 128.3, "I1" : 197.5 }, 
    "steve" : {'H1' : 160.0, 'L1' : 160.0, 'V1': 160.0, "I1" : 160.0 }, 
    "ET" : {'H1' : 2761.4017, 'L1' : 2761.4017, 'ET': 2678.2439 }
}

# sigma_t at SNR=10 in LIGO.  So, it's just 1/(20 pi sigma_f for LIGO.  
# But 1/(20 pi sigma_f)(r_ligo/r_virgo) for Virgo; 5 seconds => no localization from det.
sigma_t_dict_all = { 
    "early" : {'H1' : 0.129/1000, 'L1' : 0.129/1000, 'V1': 0.760/1000, "I1" : 5 }, 
    "mid" :   {'H1' : 0.153/1000, 'L1' : 0.153/1000, 'V1': 0.502/1000, "I1" : 5 }, 
    "late" :  {'H1' : 0.165/1000, 'L1' : 0.165/1000, 'V1': 0.353/1000, "I1" : 5 }, 
    "final" : {'H1' : 0.136/1000, 'L1' : 0.136/1000, 'V1': 0.165/1000, "I1" : 5 }, 
    "india" : {'H1' : 0.136/1000, 'L1' : 0.136/1000, 'V1': 0.165/1000, "I1" : 0.136/1000 }, 
    "steve" : {'H1' : 0.42/1000, 'L1' : 0.42/1000, 'V1': 0.42/1000, "I1" : 0.42/1000 },
    "ET" : {'H1' : 0.27463/1000, 'L1' : 0.27463/1000, 'ET': 0.163/1000 }
}






for g in range(N):
    
    g = g
    
    print "\n\n\nRun ", g, "\n"
    
    if len(sys.argv) == 2:
        net_state = sys.argv[1]
    else:
        
        print >> sys.stderr, "Usage: " + sys.argv[0] + " [" + "|".join(range_dict_all.keys()) + "]"
        print >> sys.stderr, len(sys.argv)
        for count in range(len(sys.argv)):
            print >> sys.stderr, "   " + sys.argv[count]
        sys.exit(1)
    
    
    location, response = detectors(g,bangalore=True)
    
    if net_state == "early":
        plot_proj = True  
        plot_hist = True 
    else:
        plot_proj = True
        plot_hist = True
        
    print net_state
    if net_state == "ET":
        ifos = ["H1", "L1", "ET"]
        my_networks = [ "HLE" ]

    
    if net_state == "india" or net_state == "steve":
        ifos = ["H1", "I1", "L1", "V1"]
        my_networks = [ "HILV", "ILV", "HLV", "HIV", "HIL" ]
    
    if net_state == "final":
        ifos = ["H1", "L1", "V1"]
        my_networks = [ "HLV" ]
        
    
    all_net = list(choices(ifos,2)) + \
        list(choices(ifos,3)) + \
        list(choices(ifos,4)) + \
        list(choices(ifos,5)) + \
        list(choices(ifos,6)) 
    
    four_or_more = list(choices(ifos,4)) + \
        list(choices(ifos,5)) + \
        list(choices(ifos,6)) 
    
    nets = []
    network = {}
    for dets in all_net:
        if ("H1" in dets) or (not "H2" in dets):
           if verbose: print "Keeping network:" + str(dets)
           net = ""
           for dett in dets:
                net = net + dett[0]
           nets.append(net)
           network[net] = dets
          
        else:
            if verbose: print "Not Keeping network:" + str(dets)
    
    
    ################################
    # set parameters based on network state
    ################################
    
    sigma_t = sigma_t_dict_all[net_state]
    bandwidth_dict = bandwidth_dict_all[net_state]
    range_dict = range_dict_all[net_state]
    
    ################################
    # open output file and record details
    ################################
    f = open("%s_localization_output.txt" % net_state,"w")
    
    f.write("%s configuration\n" % net_state )
    f.write("----------------------\n" )
    f.write("\nDetector ranges:\n")
    for ifo in ifos:
        f.write("%s: %d Mpc\n" % (ifo, range_dict[ifo]))
    f.write("\nSignal bandwidth:\n")
    for ifo in ifos: 
        f.write("%s: %.1f Hz\n" % (ifo, bandwidth_dict[ifo]))
    
    f.write("\nTiming accuracy for SNR 10 signal in LIGO:\n")
    for ifo in ifos:
        f.write("%s: %.2f ms\n" % (ifo, (1000 * sigma_t[ifo]) ))
    net_compare = [ my_networks ]
    
    
    ################################
    ## Time delays between sites
    ################################
    if verbose:
        print "Time delays between detectors:"
        for k1 in sort(location.keys()):
            for k2 in sort(location.keys()):
                if k1 < k2 and k1 != "H2":
                    print("IFOS: %s, %s; Light travel time %.2f ms" % (k1, k2, 
                           linalg.norm(location[k1] - location[k2])/3e8 * 1000) ) 
        print
    
    ################################
    ## Localization of Networks
    ################################
    
    
    evec = {}
    sigma = {}
    M = {}
    
    for dets in nets:
        M[dets] = localization(network[dets], sigma_t)
        evec[dets], sigma[dets] = evec_sigma(M[dets])
    
    if print_net:
        for det in sigma_t.keys():
            print("Timing accuracy in %s detector = %.2f ms\n" % 
              (det, 1000 * sigma_t[det]) )
        for net in nets:
            if len(net) > 2:
                f.write("\nFor " + str(net) + " network\n")
                f.write("Localization in degrees in eigen-directions\n")
                for i in xrange(3):
                    f.write("Direction (%.2f, %.2f, %.2f); %.2f degrees)\n" %
                (evec[net][i,0], evec[net][i,1], evec[net][i,2], 
                180/pi * sigma[net][i]) )
                a_best, a_worst = areas(M[net])
                f.write("Best case area 90: %.2f deg^2\n" % a_best)
                f.write("Worst case area 90: %.2f deg^2\n" % a_worst)
                f.write("\n")
    
    
    ################################################################
    ## Run a set of fake events, fixed distance, face on
    ################################################################
    # set up trials
    dra, ddec = pi/8.0, pi/16.0
    [ra,dec] = mgrid[-pi+dra:pi:dra,-pi/2 + ddec:pi/2:ddec]
    ra = ra.flatten()
    ra = concatenate( (ra, [0, 0]) )
    dec = dec.flatten()
    dec = concatenate( (dec, [-pi/2, pi/2]) )
    psi = 0
    cosi = 1
    loc =  xyz(ra, dec)
    
    # calculate max range in network
    rmax = max(range_dict.values())
    
    """
    if plot_proj:
      for D in [rmax, 80, 100, 160]:
          for net in my_networks:
              pylab.figure()
              pylab.axes(projection="hammer")
              pylab.grid()
              
              dets = network[net]
              for trial in xrange(len(ra)):
                  source = loc[:,trial] 
                  M, snrsq, num_found = \
                    loc_source2(dets, D, ra[trial], dec[trial], 0, cosi, response, bandwidth_dict,
                     range_dict)
          if (snrsq > 12.**2) and (num_found >= 2): plot_projection(M, source)
          else: pylab.plot(degrees(ra[trial]),degrees(dec[trial]),'rx')
          pylab.title(net)
          pylab.savefig("%s_%.0d_Mpc_localization_%s.png" %
               (net_state, D, net),  pad_inches = 0)
    """
    ################################################################
    ## Run a set of fake events, all sky all orientations
    ################################################################
    
    # Simulate a set of signals
    # Need to simulate: D, ra, dec, psi, iota
    
    D_h = 2.26 * rmax 
    
    D_max = D_h * 8 / (12/sqrt(len(ifos))) # nothing outside this can possibly give
    # combined snr > 12.
    D = random.uniform(0,1,ntrials)**(1./3) * D_max
    ra = random.uniform(0, 2*math.pi,ntrials)
    dec = arcsin( random.uniform(-1, 1,ntrials) )
    loc =  xyz(ra, dec)
    
    # calculate volume sampled and normalization factor
    vol = 4 * pi / 3 * D_max**3
    rate = 1e-6
    annual_num = vol * rate
    scaling = annual_num / ntrials
    
    psi = random.uniform(0, 2 * math.pi, ntrials)
    cosi = random.uniform(-1, 1, ntrials) 
    
    found = {}
    for ifo in ifos:
        found[ifo] = 0
    for trial in xrange(ntrials):
        for ifo in ifos:
            sig, snr = sigmasq_t(ifo, D[trial], ra[trial], dec[trial], \
                psi[trial], cosi[trial], response)
            if snr > 8.0: found[ifo] += 1
    
    tot = 0.
    for ifo in ifos:
        tot += found[ifo]
    tot /= len(ifos)
    
    f.write("----------------------------\n" )
    f.write("Ran a Monte Carlo of %g uniformly distributed sources within %.1f Mpc\n"
        % (ntrials, D_max) )
    f.write("Using a source rate of %.2g per Mpc^3 yr,\n"
            "We expect %.2f signals per year in this volume\n" 
        % (rate, annual_num) )
    f.write("So, each trial counts as %.2g signals\n" % scaling)
    
    area = {}
    for net in my_networks: area[net] = []
    
    worst = 0
    
    if plot_hist:
        for net in my_networks:
            dets = network[net]
            print >> sys.stdout, "\nThe network used is ", net
            for trial in xrange(ntrials):
                M, snrsq, num_found = loc_source2(dets, D[trial], ra[trial], dec[trial], \
                    psi[trial], cosi[trial], response, bandwidth_dict, range_dict)
                PMP = project(M, xyz(ra[trial], dec[trial])) 

                evec, sigma = evec_sigma(PMP)

                if (snrsq > 12.**2) and (num_found >= 2):
                    # event is found
                    a = 2*math.pi*2.3* (180/math.pi)**2 * sigma[0] * sigma[1]
                    area[net].append(a)
            area[net] = array(area[net])
     
            f.write("\n") 
            f.write( "\n%s network\n" % net)
            f.write("--------------\n" )
            f.write("Found %.2f signals\n" % (len(area[net]) * scaling) )
            f.write("Localized %.2f signals to within 1 deg2\n" %  
            (sum(area[net] < 1)  * scaling ) )
            f.write("Localized %.2f signals to within 2 deg2\n" %  
            (sum(area[net] < 2)  * scaling ) )
            f.write("Localized %.2f signals to within 3 deg2\n" %  
            (sum(area[net] < 3)  * scaling ) )
            f.write("Localized %.2f signals to within 5 deg2\n" %  
            (sum(area[net] < 5)  * scaling ) )
            f.write("Localized %.2f signals to within 7 deg2\n" %  
            (sum(area[net] < 7)  * scaling ) )
            f.write("Localized %.2f signals to within 10 deg2\n" %  
            (sum(area[net] < 10)  * scaling ) )
            f.write("Localized %.2f signals to within 14 deg2\n" %  
            (sum(area[net] < 14)  * scaling ) )
            f.write("Localized %.2f signals to within 16 deg2\n" %  
            (sum(area[net] < 16)  * scaling ) )
            f.write("Localized %.2f signals to within 20 deg2\n" %  
            (sum(area[net] < 20)  * scaling ) )
            f.write("Localized %.2f signals to within 24 deg2\n" %  
            (sum(area[net] < 24)  * scaling ) )
            f.write("Localized %.2f signals to within 30 deg2\n" %  
            (sum(area[net] < 30)  * scaling ) )
            f.write("Localized %.2f signals to within 39 deg2\n" %  
            (sum(area[net] < 39)  * scaling ) )
            f.write("Localized %.2f signals to within 50 deg2\n" %  
            (sum(area[net] < 50)  * scaling ) )
            f.write("Localized %.2f signals to within 70 deg2\n" %  
            (sum(area[net] < 70)  * scaling ) )
            f.write("Localized %.2f signals to within 100 deg2\n" %  
            (sum(area[net] < 100)  * scaling ) )
            f.write("Localized %.2f signals to within 140 deg2\n" %  
            (sum(area[net] < 140)  * scaling ) )
            f.write("Localized %.2f signals to within 200 deg2\n" %  
            (sum(area[net] < 200)  * scaling ) )
            f.write("Localized %.2f signals to within 300 deg2\n" %  
            (sum(area[net] < 300)  * scaling ) )
            f.write("Localized %.2f signals to within 500 deg2\n" %  
            (sum(area[net] < 500)  * scaling ) )
            f.write("Localized %.2f to worse than 500 deg2\n" % 
            (sum(area[net] > 500) * scaling)) 
            f.write("\nBest localized signal %.2f deg^2\n" % min(area[net]))
            f.write("Median localization %.2f deg^2\n" % median(area[net]) )
            f.write("Worst localized signal %.2f deg^2\n" % max(area[net]))
            # set the maximum area to 51 square degrees as we cut off the plots at 50
            area[net][area[net]>50] = 51
        print "Latitude: ", latitude[g]
        print "Longitude: ", longitude[g]
        print "Localized %.2f signals to within 5 deg2\n",  (sum(area[net] < 5)  * scaling)
       
    
    h = open("HLE localisation results 20 deg n=30 ntrials=1e5.txt", "a")    
    h.write("%.2f\t" %  (longitude[g]) )    
    h.write("%.2f\t" %  (latitude[g]) )
    h.write("%.2f\n" %  
                (sum(area[net] < 20)  * scaling ) )
    h.close()
