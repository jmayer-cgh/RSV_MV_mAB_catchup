## Core equations for transitions between compartments:
# spring birth cohort
update(M_sp) <- M_sp - n_MS_sp
update(S_sp) <- S_sp + n_MS_sp - n_SR_sp
update(R_sp) <- R_sp + n_SR_sp

# summer birth cohort
update(M_sm) <- M_sm - n_MS_sm
update(S_sm) <- S_sm + n_MS_sm - n_SR_sm
update(R_sm) <- R_sm + n_SR_sm

# autumn birth cohort
update(M_au) <- M_au - n_MS_au
update(S_au) <- S_au + n_MS_au - n_SR_au
update(R_au) <- R_au + n_SR_au

# winter birth cohort
update(M_wt) <- M_wt - n_MS_wt
update(S_wt) <- S_wt + n_MS_wt - n_SR_wt
update(R_wt) <- R_wt + n_SR_wt

# Total seroprevalence
update(R_all) <- 0.26*R_sp + 
  0.29*R_sm + 
  0.24*R_au + 
  0.20*R_wt

## Individual probabilities of transition:
# spring birth cohort
n_MS_sp <- mu * dt * M_sp # M to S
n_SR_sp <- lambda_sp * dt * S_sp # S to R

# summer birth cohort
n_MS_sm <- mu * dt * M_sm # M to S
n_SR_sm <- lambda_sm * dt * S_sm # S to R

# autumn birth cohort
n_MS_au <- mu * dt * M_au # M to S
n_SR_au <- lambda_au * dt * S_au # S to R

# winter birth cohort
n_MS_wt <- mu * dt * M_wt # M to S
n_SR_wt <- lambda_wt * dt * S_wt# S to R

## Building the FOI
# Define booleans to know which season the cohort is in
spring_FOI_sp <- if ((t <= 30.41*1.5) ||  
                     ((t > 30.41*10.5) && (t <= (30.41*13.5))) ||
                     ((t > 30.41*22.5) && (t <= (30.41*25.5))) ||
                     ((t > 30.41*34.5) && (t <= (30.41*37.5))) || 
                     ((t > 30.41*46.5) && (t <= (30.41*49.5))) || 
                     ((t > 30.41*58.5) && (t <= (30.41*61.5)))) 1 else 0
summer_FOI_sp <- if ((t > 30.41*1.5 && t <= 30.41*4.5 ) ||                      # FOI in summer for those born in spring
                     ((t > 30.41*13.5) && (t <= 30.41*16.5)) ||
                     ((t > 30.41*25.5) && (t <= 30.41*28.5)) ||
                     ((t > 30.41*37.5) && (t <= 30.41*40.5)) ||
                     ((t > 30.41*49.5) && (t <= 30.41*52.5))) 1 else 0 
autumn_FOI_sp <- if ( (t > 30.41*4.5 && t <= 30.41*7.5 ) ||                    
                      ((t > 30.41*16.5) && (t <= 30.41*19.5)) ||
                      ((t > 30.41*28.5) && (t <= 30.41*31.5)) ||
                      ((t > 30.41*40.5) && (t <= 30.41*43.5)) ||
                      ((t > 30.41*52.5) && (t <= 30.41*55.5))) 1 else 0
winter_FOI_sp <- if ( (t > 30.41*7.5 && t <= 30.41*10.5 ) ||                   
                      ((t > 30.41*19.5) && (t <= 30.41*22.5)) ||
                      ((t > 30.41*31.5) && (t <= 30.41*34.5)) ||
                      ((t > 30.41*43.5) && (t <= 30.41*46.5)) ||
                      ((t > 30.41*55.5) && (t <= 30.41*58.5))) 1 else 0
spring_FOI_sm <- if ( (t > 30.41*7.5 && t <= 30.41*10.5 ) ||                   
                      ((t > 30.41*19.5) && (t <= 30.41*22.5)) ||
                      ((t > 30.41*31.5) && (t <= 30.41*34.5)) ||
                      ((t > 30.41*43.5) && (t <= 30.41*46.5)) ||
                      ((t > 30.41*55.5) && (t <= 30.41*58.5))) 1 else 0
summer_FOI_sm <- if ((t <= 30.41*1.5) ||  
                     ((t > 30.41*10.5) && (t <= (30.41*13.5))) ||
                     ((t > 30.41*22.5) && (t <= (30.41*25.5))) ||
                     ((t > 30.41*34.5) && (t <= (30.41*37.5))) || 
                     ((t > 30.41*46.5) && (t <= (30.41*49.5))) || 
                     ((t > 30.41*58.5) && (t <= (30.41*61.5)))) 1 else 0
autumn_FOI_sm <- if ((t > 30.41*1.5 && t <= 30.41*4.5 ) ||                      # FOI in autumn for those born in summer
                     ((t > 30.41*13.5) && (t <= 30.41*16.5)) ||
                     ((t > 30.41*25.5) && (t <= 30.41*28.5)) ||
                     ((t > 30.41*37.5) && (t <= 30.41*40.5)) ||
                     ((t > 30.41*49.5) && (t <= 30.41*52.5))) 1 else 0
winter_FOI_sm <-  if ( (t > 30.41*4.5 && t <= 30.41*7.5 ) ||                    
                       ((t > 30.41*16.5) && (t <= 30.41*19.5)) ||
                       ((t > 30.41*28.5) && (t <= 30.41*31.5)) ||
                       ((t > 30.41*40.5) && (t <= 30.41*43.5)) ||
                       ((t > 30.41*52.5) && (t <= 30.41*55.5))) 1 else 0 
spring_FOI_au <-  if ( (t > 30.41*4.5 && t <= 30.41*7.5 ) ||                    
                       ((t > 30.41*16.5) && (t <= 30.41*19.5)) ||
                       ((t > 30.41*28.5) && (t <= 30.41*31.5)) ||
                       ((t > 30.41*40.5) && (t <= 30.41*43.5)) ||
                       ((t > 30.41*52.5) && (t <= 30.41*55.5))) 1 else 0
summer_FOI_au <- if ( (t > 30.41*7.5 && t <= 30.41*10.5 ) ||                   
                      ((t > 30.41*19.5) && (t <= 30.41*22.5)) ||
                      ((t > 30.41*31.5) && (t <= 30.41*34.5)) ||
                      ((t > 30.41*43.5) && (t <= 30.41*46.5)) ||
                      ((t > 30.41*55.5) && (t <= 30.41*58.5))) 1 else 0
autumn_FOI_au <- if ((t <= 30.41*1.5) ||  
                     ((t > 30.41*10.5) && (t <= (30.41*13.5))) ||
                     ((t > 30.41*22.5) && (t <= (30.41*25.5))) ||
                     ((t > 30.41*34.5) && (t <= (30.41*37.5))) || 
                     ((t > 30.41*46.5) && (t <= (30.41*49.5))) || 
                     ((t > 30.41*58.5) && (t <= (30.41*61.5)))) 1 else 0
winter_FOI_au <- if ((t > 30.41*1.5 && t <= 30.41*4.5 ) ||                      
                     ((t > 30.41*13.5) && (t <= 30.41*16.5)) ||
                     ((t > 30.41*25.5) && (t <= 30.41*28.5)) ||
                     ((t > 30.41*37.5) && (t <= 30.41*40.5)) ||
                     ((t > 30.41*49.5) && (t <= 30.41*52.5))) 1 else 0
spring_FOI_wt <- if ((t > 30.41*1.5 && t <= 30.41*4.5 ) ||                      
                     ((t > 30.41*13.5) && (t <= 30.41*16.5)) ||
                     ((t > 30.41*25.5) && (t <= 30.41*28.5)) ||
                     ((t > 30.41*37.5) && (t <= 30.41*40.5)) ||
                     ((t > 30.41*49.5) && (t <= 30.41*52.5))) 1 else 0
summer_FOI_wt <-  if ( (t > 30.41*4.5 && t <= 30.41*7.5 ) ||                    
                       ((t > 30.41*16.5) && (t <= 30.41*19.5)) ||
                       ((t > 30.41*28.5) && (t <= 30.41*31.5)) ||
                       ((t > 30.41*40.5) && (t <= 30.41*43.5)) ||
                       ((t > 30.41*52.5) && (t <= 30.41*55.5))) 1 else 0
autumn_FOI_wt <- if ( (t > 30.41*7.5 && t <= 30.41*10.5 ) ||                   
                      ((t > 30.41*19.5) && (t <= 30.41*22.5)) ||
                      ((t > 30.41*31.5) && (t <= 30.41*34.5)) ||
                      ((t > 30.41*43.5) && (t <= 30.41*46.5)) ||
                      ((t > 30.41*55.5) && (t <= 30.41*58.5))) 1 else 0
winter_FOI_wt <- if ((t <= 30.41*1.5) ||  
                     ((t > 30.41*10.5) && (t <= (30.41*13.5))) ||
                     ((t > 30.41*22.5) && (t <= (30.41*25.5))) ||
                     ((t > 30.41*34.5) && (t <= (30.41*37.5))) || 
                     ((t > 30.41*46.5) && (t <= (30.41*49.5))) || 
                     ((t > 30.41*58.5) && (t <= (30.41*61.5)))) 1 else 0

# Putting it all together into four FOIs
lambda_sp = (summer_comp + spring_comp) * spring_FOI_sp + 
  summer_comp * summer_FOI_sp + 
  (summer_comp + autumn_comp) * autumn_FOI_sp +
  (summer_comp + winter_comp) * winter_FOI_sp 

lambda_sm = (summer_comp + spring_comp) * spring_FOI_sm + 
  summer_comp * summer_FOI_sm + 
  (summer_comp + autumn_comp) * autumn_FOI_sm + 
  (summer_comp + winter_comp) * winter_FOI_sm 

lambda_au = (summer_comp + spring_comp) * spring_FOI_au + 
  summer_comp * summer_FOI_au + 
  (summer_comp + autumn_comp) * autumn_FOI_au + 
  (summer_comp + winter_comp) * winter_FOI_au 

lambda_wt = (summer_comp + spring_comp) * spring_FOI_wt + 
  summer_comp * summer_FOI_wt +
  (summer_comp + autumn_comp) * autumn_FOI_wt + 
  (summer_comp + winter_comp) * winter_FOI_wt 

## Initial states:
# spring cohort
initial(M_sp) <- prop * M_sp_ini # only a proportion of children is born protected
initial(S_sp) <- (1-prop) * M_sp_ini
initial(R_sp) <- R_sp_ini

# summer cohort
initial(M_sm) <- prop * M_sm_ini
initial(S_sm) <- (1-prop) * M_sm_ini
initial(R_sm) <- R_sm_ini

# autumn cohort
initial(M_au) <- prop * M_au_ini
initial(S_au) <- (1-prop) * M_au_ini
initial(R_au) <- R_au_ini

# winter cohort
initial(M_wt) <- prop * M_wt_ini
initial(S_wt) <- (1-prop) * M_wt_ini
initial(R_wt) <- R_wt_ini

# Total
initial(R_all) <- 0.26*R_sp_ini + 
  0.29*R_sm_ini + 
  0.24*R_au_ini + 
  0.20*R_wt_ini

## User defined parameters - default in parentheses:
M_sp_ini <- parameter(1 - 1e-12)
R_sp_ini <- parameter(1e-12)

M_sm_ini <- parameter(1 - 1e-12)
R_sm_ini <- parameter(1e-12)

M_au_ini <- parameter(1 - 1e-12)
R_au_ini <- parameter(1e-12)

M_wt_ini <- parameter(1 - 1e-12)
R_wt_ini <- parameter(1e-12)


# transition parameters
spring_comp <- parameter(1e-05)
summer_comp <- parameter(0.02002)
autumn_comp <- parameter(3e-05)
winter_comp <- parameter(4e-05)
mu <- parameter (0.0050) # parameter(1/59.50121)

# Proportion born with maternal immunity
prop <- parameter(1)

# Comparison function
# By season
N_spring <- data()
n_infection_spring <- data()
n_infection_spring ~ Binomial(N_spring, R_sp)

N_summer <- data()
n_infection_summer <- data()
n_infection_summer ~ Binomial(N_summer, R_sm)

N_autumn <- data()
n_infection_autumn <- data()
n_infection_autumn ~ Binomial(N_autumn, R_au)

N_winter <- data()
n_infection_winter <- data()
n_infection_winter ~ Binomial(N_winter, R_wt)