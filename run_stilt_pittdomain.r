#!/usr/bin/env Rscript
# STILT R Executable
# For documentation, see https://uataq.github.io/stilt/
# Ben Fasoli

library("optparse")
 
option_list = list(
  make_option(c("-z", "--alt"), type="double", default=NULL, 
              help="Receptor Altitude (m)", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Receptor Filename (.csv)", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$alt) + is.null(opt$file)){
  print_help(opt_parser)
  stop("The receptor altitude and file name must be supplied.n", call.=FALSE)
}

# User inputs ------------------------------------------------------------------
project <- 'halo'
stilt_wd <- file.path('/discover/nobackup/smcrowe1/stilt', project)
output_wd <- file.path(stilt_wd,paste0(substr(opt$file,11,nchar(opt$file)-4),'_',format(opt$alt,nsmall=1)))
lib.loc <- .libPaths()[1]

print(paste0("Output Dir:", output_wd))

# Parallel simulation settings
n_cores <- 36 
n_nodes <- 1
slurm   <- n_nodes > 1
slurm_options <- list(
  time      = '08:00:00',
  account   = 'smcrowe1',
  partition = 'compute'
)

# Simulation timing, yyyy-mm-dd HH:MM:SS (UTC)
##t_start <- '2015-12-10 00:00:00'
##t_end   <- '2015-12-10 00:00:00'
##run_times <- seq(from = as.POSIXct(t_start, tz = 'UTC'),
##                 to   = as.POSIXct(t_end, tz = 'UTC'),
##                 by   = 'hour')

# Receptor location(s)
##lati <- 40.5
##long <- -112.0
##zagl <- 5

# Expand the run times, latitudes, and longitudes to form the unique receptors
# that are used for each simulation
##receptors <- expand.grid(run_time = run_times, lati = lati, long = long,
##                         zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)

receptors <- read.csv(opt$file)
receptors$run_time <- as.POSIXct(paste0(receptors$UTC_date, receptors$UTC_time), tz = "UTC")
receptors$zagl = opt$alt
n_receptors = nrow(receptors)

#to_add = receptors[1:n_receptors,]
#to_add$long = to_add$long - (i*0.00593469134*4) #2km
#to_add$long = as.numeric(sprintf("%.5f",to_add$long))

#receptors = rbind(receptors, to_add)

browser()

# Footprint grid settings, must set at least xmn, xmx, ymn, ymx below
hnf_plume <- T
projection <- '+proj=longlat'
smooth_factor <- 1
time_integrate <- F
xmn <- -75.7
xmx <- -72.1
ymn <- 39.2
ymx <- 42.0
xres <- 0.01
yres <- 0.01

# Meteorological data input
met_path           <- '/discover/nobackup/smcrowe1/hrrr/'
met_file_format    <- '%Y%m%d_%Hz_hrrr'
met_subgrid_buffer <- 0.2
met_subgrid_enable <- F
met_subgrid_levels <- NA
n_met_min          <- 1

# Model control
n_hours       <- -24
numpar        <- 1000
rm_dat        <- T
run_foot      <- T
run_trajec    <- T
simulation_id <- NA
timeout       <- 3600
varsiwant  <- c('time', 'indx', 'long', 'lati', 'zagl', 'sigw', 'tlgr', 'zsfc',
                'icdx', 'temp', 'samt', 'foot', 'shtf', 'tcld', 'dmas', 'dens',
                'rhfr', 'sphu', 'lcld', 'zloc', 'dswf', 'wout', 'mlht', 'rain',
                'crai', 'pres', 'whtf', 'temz', 'zcfx')

# Transport and dispersion settings
capemin     <- -1
cmass       <- 0
conage      <- 48
cpack       <- 1
delt        <- 1
dxf         <- 1
dyf         <- 1
dzf         <- 0.01
efile       <- ''
emisshrs    <- 0.01
frhmax      <- 3
frhs        <- 1
frme        <- 0.1
frmr        <- 0
frts        <- 0.1
frvs        <- 0.01
hscale      <- 10800
ichem       <- 8
idsp        <- 2
initd       <- 0
k10m        <- 1
kagl        <- 1
kbls        <- 1
kblt        <- 5
kdef        <- 0
khinp       <- 0
khmax       <- 9999
kmix0       <- 150
kmixd       <- 3
kmsl        <- 0
kpuff       <- 0
krand       <- 4
krnd        <- 6
kspl        <- 1
kwet        <- 1
kzmix       <- 0
maxdim      <- 1
maxpar      <- numpar
mgmin       <- 10
mhrs        <- 9999
nbptyp      <- 1
ncycl       <- 0
ndump       <- 0
ninit       <- 1
nstr        <- 0
nturb       <- 0
nver        <- 0
outdt       <- 0
p10f        <- 1
pinbc       <- ''
pinpf       <- ''
poutf       <- ''
qcycle      <- 0
rhb         <- 80
rht         <- 60
splitf      <- 1
tkerd       <- 0.18
tkern       <- 0.18
tlfrac      <- 0.1
tout        <- 0
tratio      <- 0.75
tvmix       <- 1
veght       <- 0.5
vscale      <- 200
vscaleu     <- 200
vscales     <- -1
wbbh        <- 0
wbwf        <- 0
wbwr        <- 0
wvert       <- FALSE
w_option    <- 0
zicontroltf <- 0
ziscale     <- rep(list(rep(1, 24)), nrow(receptors))
z_top       <- 25000

# Transport error settings
horcoruverr <- NA
siguverr    <- NA
tluverr     <- NA
zcoruverr   <- NA

horcorzierr <- NA
sigzierr    <- NA
tlzierr     <- NA


# Interface to mutate the output object with user defined functions
before_trajec <- function() {output}
before_footprint <- function() {output}


# Source dependencies ----------------------------------------------------------
setwd(stilt_wd)
source('r/dependencies.r')


# Structure out directory ------------------------------------------------------
# Outputs are organized in three formats. by-id contains simulation files by
# unique simulation identifier. particles and footprints contain symbolic links
# to the particle trajectory and footprint files in by-id
system(paste0('rm -r ', output_wd, '/footprints'), ignore.stderr = T)
if (run_trajec) {
  system(paste0('rm -r ', output_wd, '/by-id'), ignore.stderr = T)
  system(paste0('rm -r ', output_wd, '/met'), ignore.stderr = T)
  system(paste0('rm -r ', output_wd, '/particles'), ignore.stderr = T)
}
for (d in c('by-id', 'particles', 'footprints')) {
  d <- file.path(output_wd, d)
  if (!file.exists(d))
    dir.create(d, recursive = T)
}


# Run trajectory simulations ---------------------------------------------------
stilt_apply(FUN = simulation_step,
            simulation_id = simulation_id,
            slurm = slurm, 
            slurm_options = slurm_options,
            n_cores = n_cores,
            n_nodes = n_nodes,
            before_footprint = list(before_footprint),
            before_trajec = list(before_trajec),
            lib.loc = lib.loc,
            capemin = capemin,
            cmass = cmass,
            conage = conage,
            cpack = cpack,
            delt = delt,
            dxf = dxf,
            dyf = dyf,
            dzf = dzf,
            efile = efile,
            emisshrs = emisshrs,
            frhmax = frhmax,
            frhs = frhs,
            frme = frme,
            frmr = frmr,
            frts = frts,
            frvs = frvs,
            hnf_plume = hnf_plume,
            horcoruverr = horcoruverr,
            horcorzierr = horcorzierr,
            hscale = hscale,
            ichem = ichem,
            idsp = idsp,
            initd = initd,
            k10m = k10m,
            kagl = kagl,
            kbls = kbls,
            kblt = kblt,
            kdef = kdef,
            khinp = khinp,
            khmax = khmax,
            kmix0 = kmix0,
            kmixd = kmixd,
            kmsl = kmsl,
            kpuff = kpuff,
            krand = krand,
            krnd = krnd,
            kspl = kspl,
            kwet = kwet,
            kzmix = kzmix,
            maxdim = maxdim,
            maxpar = maxpar,
            met_file_format = met_file_format,
            met_path = met_path,
            met_subgrid_buffer = met_subgrid_buffer,
            met_subgrid_enable = met_subgrid_enable,
            met_subgrid_levels = met_subgrid_levels,
            mgmin = mgmin,
            n_hours = n_hours,
            n_met_min = n_met_min,
            ncycl = ncycl,
            ndump = ndump,
            ninit = ninit,
            nstr = nstr,
            nturb = nturb,
            numpar = numpar,
            nver = nver,
            outdt = outdt,
            output_wd = output_wd,
            p10f = p10f,
            pinbc = pinbc,
            pinpf = pinpf,
            poutf = poutf,
            projection = projection,
            qcycle = qcycle,
            r_run_time = receptors$run_time,
            r_lati = receptors$lati,
            r_long = receptors$long,
            r_zagl = receptors$zagl,
            rhb = rhb,
            rht = rht,
            rm_dat = rm_dat,
            run_foot = run_foot,
            run_trajec = run_trajec,
            siguverr = siguverr,
            sigzierr = sigzierr,
            smooth_factor = smooth_factor,
            splitf = splitf,
            stilt_wd = stilt_wd,
            time_integrate = time_integrate,
            timeout = timeout,
            tkerd = tkerd,
            tkern = tkern,
            tlfrac = tlfrac,
            tluverr = tluverr,
            tlzierr = tlzierr,
            tout = tout,
            tratio = tratio,
            tvmix = tvmix,
            varsiwant = list(varsiwant),
            veght = veght,
            vscale = vscale,
            vscaleu = vscaleu,
            vscales = vscales,
            w_option = w_option,
            wbbh = wbbh,
            wbwf = wbwf,
            wbwr = wbwr,
            wvert = wvert,
            xmn = xmn,
            xmx = xmx,
            xres = xres,
            ymn = ymn,
            ymx = ymx,
            yres = yres,
            zicontroltf = zicontroltf,
            ziscale = ziscale,
            z_top = z_top,
            zcoruverr = zcoruverr)
