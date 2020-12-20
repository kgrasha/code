import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

comp = fits.open('/Users/grasha/Research/typhoon/data/N1365_1_comp.fits')

ha = comp['HALPHA'].data[0,:,:]
haerr = comp['HALPHA_ERR'].data[0,:,:]
hb = comp['HBETA'].data[0,:,:]
hberr = comp['HBETA_ERR'].data[0,:,:]
o3 = comp['OIII5007'].data[0,:,:]
o3err = comp['OIII5007_ERR'].data[0,:,:]
n2 = comp['NII6583'].data[0,:,:]
n2err = comp['NII6583_ERR'].data[0,:,:]

sii16 = comp['SII6716'].data[0,:,:]
sii31 = comp['SII6731'].data[0,:,:]
s16err = comp['SII6716_ERR'].data[0,:,:]
s31err = comp['SII6731_ERR'].data[0,:,:]
oii26 = comp['OII3726'].data[0,:,:]
oii29 = comp['OII3729'].data[0,:,:]
oii26err = comp['OII3726_ERR'].data[0,:,:]
oii29err = comp['OII3729_ERR'].data[0,:,:]

# add together the SII and OII lines
s2 = sii16 + sii31
s2err = np.sqrt(s16err**2 + s31err**2)

o2 = oii26 + oii29
o2err = np.sqrt(oii26err**2 + oii29err**2)

# set snr ratio
snr = 3.0

# thow out spaxels with S/N < 3 for Hbeta, SII, N2, and Halpha
for i in range(0, ha.shape[0]):
    for j in range(0, ha.shape[1]):
        if ha[i][j]/haerr[i][j] <=snr or hb[i][j]/hberr[i][j] <= snr or s2[i][j]/s2err[i][j] <= snr or n2[i][j]/n2err[i][j] <= snr:
            ha[i][j] = np.nan
            hb[i][j] = np.nan
            o3[i][j] = np.nan
            n2[i][j] = np.nan
            s2[i][j] = np.nan
            o2[i][j] = np.nan


# thow out spaxels with flux values<0 for Hbeta, SII, N2, Halpha, OIII
for i in range(0, ha.shape[0]):
    for j in range(0, ha.shape[1]):
        if ha[i][j]<=0 or hb[i][j] <= 0 or s2[i][j]<=0 or n2[i][j]== 0 or o3[i][j]<=0:
            ha[i][j] = np.nan
            hb[i][j] = np.nan
            o3[i][j] = np.nan
            n2[i][j] = np.nan
            s2[i][j] = np.nan
            o2[i][j] = np.nan


# N2S2-N2Ha Dopita+2016 metallicty diagnostic
y = np.log10(n2/s2) + 0.264*np.log10(n2/ha)
met_d16 = 8.77 + y + 0.45*(y+0.3)**5

# Throw out spaxesl that have initital guesses at a range we don't want
for i in range(0, met_d16.shape[0]):
    for j in range(0, met_d16.shape[1]):
        if met_d16[i][j] > 9.4 or met_d16[i][j] < 7.5:
            met_d16[i][j] = np.nan


# look at the image, make sure the range looks good
plt.imshow(met_d16, origin='lower', cmap=plt.get_cmap('plasma'))
plt.ylim(100,500)
plt.title(r'N2S2 - N2H$\alpha$')
ax = plt.colorbar()
ax.set_label(r'12 + log(O/H)',fontsize=20)
plt.show()

# look at the histogram, make sure the range looks good
plt.hist(met_d16.flatten(), bins=30, edgecolor='black')
plt.xlabel('12+log(O/H)')
plt.ylabel('Number')
plt.show()



# Use the N2S2 as an initial guess for logU from Kewley+19 ARAA
zz = met_d16 #make an initial guess for metallicity


# Need to iterate to get U and metallicity from Kewley+2019 ARAA
for i in range(6):
    print(i)

    # set up eqns to calculate U
    A = 13.768
    B = 9.4940
    C = -4.3223
    D = -2.3531
    E = -0.5769
    F = 0.2794
    G = 0.1574
    H = 0.0890
    I = 0.0311
    J = 0.0

    x = np.log10(o3/o2)

    y = zz 
    UO3O2 = A + B*x + C*y + D*x*y + E*x**2 + F*y**2 + G*x*y**2 + H*y*x**2 + I*x**3 + J*y**3
    print(UO3O2[366,49]) #print U at a random pixel

    #set up equations to calculate met (N2Ha) using the U calculated above
    A = 10.526 
    B = 1.9958
    C = -0.6741
    D = 0.2892
    E = 0.5712
    F = -0.6597
    G = 0.0101
    H = 0.0800
    I = 0.0782
    J = -0.0982

    x = np.log10(n2/ha)

    y = UO3O2
    zN2Ha = A + B*x + C*y + D*x*y + E*x**2 + F*y**2 + G*x*y**2 + H*y*x**2 + I*x**3 + J*y**3
    zz = zN2Ha
    print(zN2Ha[366,49], 'metallicity at iteration', i, '\n') #print met at the random pixel

#compare the printed metallicity values at for the iterations. Have the values stablized to ~0.01 dex between the last two iterations?

# You'll notice that the center of the galaxy isn't constrained very well and have U and metallicity values way above the accepted range of U and metallicity. The main reason for this is that this region of the galaxy isn't being powered by star formation but is instead powered by shocks/AGNs. What we can do is pull up the BPT diagrams and exclude spaxels that are not consistent with being powered by star formation (i.e., they will lie over the Kauffmann star formation line in the BPT-NII plot).
# https://sites.google.com/site/agndiagnostics/home/bpt
# If you go back to the beginning of the program and plot the BPT-NII plot (OIII/Hb vs NII/Ha) and exclude spaxels over the Kauffmann+2003 line, does this improve your range of U?


plt.figure()
plt.imshow(UO3O2, origin='lower', cmap=plt.get_cmap('plasma'))
plt.ylim(100,500)
ax = plt.colorbar()
ax.set_label(r'log U',fontsize=20)
plt.show()

#Look at the image for the N2Ha metallicity distrubution
plt.figure()
plt.imshow(zN2Ha, origin='lower', cmap=plt.get_cmap('plasma'))
plt.ylim(100,500)
plt.title(r'N2H$\alpha$')
plt.colorbar()
plt.show()


# compare the distribution for N2S2 to N2Ha

plt.hist(zN2Ha.flatten(), alpha=0.6, edgecolor='black', color='teal', label='N2Hα', bins=30)
plt.hist(met_d16.flatten(), alpha=0.5, edgecolor='black', color='mediumvioletred', label='N2S2', bins=30)
plt.legend()
plt.xlabel('12+log(O/H)')
plt.ylabel('Number')
plt.show()

# You will notice that N2Ha histogram is very similar in shape to the N2S2 diagnostic. This is good. There is a tail for the N2Ha distribution that extends to much larger values than we want. If you compare this with the image for the N2Ha metallicity distrubution, these high values are all in the centre of the galaxy, which has an AGN source which is messing up the center values. If you use the BPT diagram and remove the spaxels in the center of the galaxy, maybe this will help improve the range of N2Ha metallicity values. 
