  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //FUNCTIONS
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
functions {
  real weibull_pdf(real y, real k, real lambda) {
    if (y < 0) return 0; // Ensure y >= 0
    return ((k / lambda) * pow(y / lambda, k - 1) * exp(-pow(y / lambda, k)));
  }

  real weibull_cdf_manual(real x, real shape, real scale) {
    if (x < 0) {
      return 0;
    } else {
      return 1 - exp(-pow(x / scale, shape));
    }
  }
  
  real poisson_pdf_manual(int x, real lambda) {
    return (exp(-lambda)*pow(lambda, x)/exp(lgamma(x + 1))); //exp(lgamma(x+1))=tgamma(x+1)
}
}
 
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //DATA
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

data {
  //-----------------------------------
  //Global Settings
  //-----------------------------------
  int<lower=1, upper=4> Adstock_Method; // Method for adstock 1: Geometric Decay Non-Weighted | 2: Geometric Decay Weighted / 3: Weibull PDF Decay / 4: Poisson Decay
  int<lower=1> M; // M: cantidad de canales (con spend)
  int<lower=0> C; // C: catidad de variales de control sin adstock ni saturation effect (Impacto directo en ventas)
  int<lower=1> C_used; // For dummie uses under C=0
  int<lower=1> T; // T: cantidad de observaciones temporales
  int<lower=0> L; // L: cantidad de lags máximos
  int<lower=0> M_sim; // M_sim: amount of simulations
  
  //-----------------------------------
  //Y - Variable Objetivo
  //-----------------------------------
  vector[T] Y ;  //Y: variable ojetivo, revenue u otro KPI
  int<lower=0,upper=1> is_test[T]; // Indicador: 1 si test, 0 si train
  
  //-----------------------------------
  //TEMP - Variables Temporales
  //-----------------------------------
  int<lower=1, upper=2> Trend_Method; // Method for trend 1: Years | 2: Months
  
  int<lower=0> t_ind[T]; 
  int<lower=0> semana_ano[T]; 
  int<lower=0> mes[T];
  int<lower=0> ano[T];
  
  //-----------------------------------
  //xi - Variables de Marketing
  //-----------------------------------
  matrix<lower=0>[T,M] x ;  //Canales con spend
  
  //-----------------------------------
  //ci - Variables Control
  //-----------------------------------
  matrix<lower=0>[T,C_used] ct ;  //Variales de control sin adstock ni saturation effect (Impacto directo en ventas)
  
  //-----------------------------------
  //HiperParámetros Prioris
  //-----------------------------------
  // Parámetros generales
  real<lower=1> beta_0_shape;
  real<lower=0> beta_0_scale;
  
  real<lower=1> sigma_e_shape;
  real<lower=0> sigma_e_scale;
  
  real tau_mean;
  real<lower=0> tau_sigma;
  
  real<lower=0> sigma_psi_mean;
  real<lower=0> sigma_psi_sigma;
  
  // Por canales (X)
  // ADSTOCK
  // Geometric Adstock
  vector<lower=0>[M] theta_x_a;
  vector<lower=0>[M] theta_x_b;
  
  // Weibull Adstock
  vector<lower=1>[M] shape_x_shape;
  vector<lower=0>[M] shape_x_scale;
  
  vector<lower=1>[M] scale_x_shape;
  vector<lower=0>[M] scale_x_scale;
  
  // Poisson Adstock
  vector<lower=0>[M] lambda_x_min;
  vector<lower=0>[M] lambda_x_max;
  
  // SATUTATION
  vector<lower=1>[M] beta_x_shape;
  vector<lower=0>[M] beta_x_scale;
  
  vector<lower=1>[M] gamma_x_shape;
  vector<lower=0>[M] gamma_x_scale;
  
  vector<lower=1>[M] alpha_x_shape;
  vector<lower=0>[M] alpha_x_scale;
  
  // Por Controles (c)
  vector[C_used] beta_ct_mu;
  vector<lower=0>[C_used] beta_ct_sigma;
  
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//PARAMETERS
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
parameters {
  
  //-----------------------------------
  //GLOBAL Variables
  //-----------------------------------
  real<lower=0> beta_0;             // Intercept
  real<lower=0> sigma_e;   // Standard deviation of the white noise around signal effect
  
  //-----------------------------------
  //TEMP - Variables Temporales
  //-----------------------------------
  real tau; // Efecto tendencia (trend) - variación de ventas debido al efecto tendencia
  real psi[12]; // Efecto estacional (seasonal) - variación de ventas debido al efecto estacional (mensual)
  
  real<lower=0> sigma_psi;  // Variabilidad en los efectos estacionales (variabilidad entre meses)
  
  //-----------------------------------
  //xi - Variables de Marketing
  //-----------------------------------
  vector<lower=0, upper=1>[M] theta_x;   // Theta for each channel ADSTOCK GEOMETRIC
  vector<lower=1>[M] shape_x;   // Shape for each channel ADSTOCK WEIBULL
  vector<lower=0>[M] scale_x;   // Scale for each channel ADSTOCK WEIBULL
  vector<lower=0>[M] lambda_x;   // Lambda for each channel ADSTOCK POISSON
  
  vector<lower=0>[M] beta_x;     // Beta for each channel
  vector<lower=0>[M] gamma_x;   // Gamma (ec) for each channel SATURATION
  
  vector<lower=0>[M] alpha_x;   // Alpha (slope) for each channel SATURATION
  
  
  //-----------------------------------
  //ci - Variables Control
  //-----------------------------------
  vector[C_used] beta_ct; //Beta for each control variable (How an increase in the Control Variable affects KPI-Y)
}

transformed parameters {
  //-----------------------------------
  // TEMPORAL EFFECTS
  //-----------------------------------
  
  //TREND
  //----------------------------------
  real trend[T];
  
  real estima_tendencia[T];
  
  if(Trend_Method==1){  //Estimación a través de los anos
      estima_tendencia = ano;
        }
        
  if(Trend_Method==2){  //Estimación a través de las semanas
      estima_tendencia = t_ind;
        }
  
  for (t in 1:T){
    trend[t] = tau*estima_tendencia[t];
  } 

  //STATIONALITY
  //----------------------------------
  real seasonal[T];
  
  for (t in 1:T){
    seasonal[t] = psi[mes[t]];
  } 

  
  //TEMPORAL OVERALL
  //----------------------------------
  real Y_temporal[T];
  
  for (t in 1:T){
    Y_temporal[t] = beta_0 + trend[t] + seasonal[t];  //Cuánto de los ingresos depende de los diferentes canales de marketing
  }
  
  //-----------------------------------
  //ADSTOCK
  //-----------------------------------

      vector[M] W_x_factor;       // Factor de conversión de resultados no ponderados a resultados ponderados
      vector<lower=0>[M] W_x_den; // Denominador del factor de conversión de resultados a ponderados
      
      matrix<lower=0>[(L+1),M] w_x_s; // Matriz de pesos para cada lag-s (en cada medio m)
      
          if(Adstock_Method == 1 || Adstock_Method == 2){  //GEOMETRICAL ADSTOCK
                  for (m in 1:M){
                  for (s in 1:(L+1)){
                    w_x_s[s,m] =  theta_x[m]^(s-1);         
                  }
                  }
                }
    
          

          if(Adstock_Method==3){  //WEIBULL PDF ADSTOCK
                  for (m in 1:M){
                    real s_mode_weibull = scale_x[m] * pow((shape_x[m] - 1) / shape_x[m], 1 / shape_x[m]); // Lag s that will return the maximum adstock impact
                    real max_weibull = weibull_pdf(s_mode_weibull, shape_x[m], scale_x[m]);
                    
                  for (s in 1:(L+1)){
                    
                    if ((shape_x[m]>1) && (s==1) ){
                        real weibull_density = 1e-10;
                        w_x_s[s,m] =  weibull_density/max_weibull; 
                      }else{
                        real weibull_density = weibull_pdf(s-1, shape_x[m], scale_x[m]);
                        w_x_s[s,m] =  weibull_density/max_weibull;}         
                  }
                  }
                }
                
                if(Adstock_Method==4){  //POISSON DISCRETE ADSTOCK
                  for (m in 1:M){
                  for (s in 1:(L+1)){
                  
                      real poisson_probability = poisson_pdf_manual(s-1, lambda_x[m]);
                      w_x_s[s,m] =  poisson_probability+1e-10;
           
                  }
                  }
                }
      
      for(m in 1:M){
          W_x_den[m]=0;                     // Valor de inicialización
              for (s in 1:(L+1)){
                W_x_den[m] +=  w_x_s[s,m];         // Para cada lag s-1, va sumando el efecto
              }
          W_x_factor[m] = 1 / W_x_den[m];
      }
      
      matrix<lower=0>[T,M] x_ads;
      matrix<lower=0>[T,M] x_pre_ads;
      vector<lower=0>[M] gamma_w_x;
      vector<lower=0>[M] gamma_non_w_x;
      
      for(m in 1:M){
          for (t in 1:T){
             
            x_pre_ads[t,m]=0;       // Valor de inicialización
              for (s in 1:(L+1)) {
                  if (t - (s-1) > 0) {
                      x_pre_ads[t,m] += w_x_s[s,m] * x[(t - (s-1)),m]; // Para cada lag s, va sumando el efecto
              }}
            
            x_ads[t,m] = W_x_factor[m]*x_pre_ads[t,m];
            }
            
            // WEIGHTED?
            if(Adstock_Method==1){  //Non-weighted Estimation
              gamma_non_w_x[m] = gamma_x[m]; 
              gamma_w_x[m]     = gamma_x[m]*W_x_factor[m];
            }
            
            if(Adstock_Method!=1){  //Weighted Estimation
              gamma_non_w_x[m] = gamma_x[m]/W_x_factor[m]; 
              gamma_w_x[m]     = gamma_x[m];
            }
       }

  //-----------------------------------
  //SATURATION - HILL
  //-----------------------------------
  matrix<lower=0>[T,M] x_sat;
  matrix<lower=0>[T,M] x_hilladstocked;
  
  for (m in 1:M){
  for (t in 1:T){
  x_sat[t,m]=1/(1+ (gamma_w_x[m]/x_ads[t,m])^alpha_x[m] );
  x_hilladstocked[t,m]= beta_x[m] * x_sat[t,m]; // Revenue € de Y[t] dependiente del medio x1 [en el momento t, debido al anuncio acumulado y saturado en ese momento t]
  }}

  //-----------------------------------
  //MARKETING EFFECTS (HILLADSTOCKED)
  //-----------------------------------
  real Y_marketing[T];
  
  for (t in 1:T){
    Y_marketing[t] = sum(x_hilladstocked[t]);  //Cuánto de los ingresos depende de los diferentes canales de marketing
  }
  
  vector[M] ROI_x;
  
  for (m in 1:M) {
    ROI_x[m]= sum(x_hilladstocked[,m])/sum(x[,m]);
  }
  
  
  //-----------------------------------
  // CONTROL EFFECTS
  //-----------------------------------
  matrix[T,C_used] ct_impact;
  real Y_control[T];
  
  for (c in 1:C_used){
  for (t in 1:T){
  ct_impact[t,c]= beta_ct[c] * ct[t,c]; // Revenue € de Y[t] dependiente del control c en el momento t
  }}
  
  for (t in 1:T){
    Y_control[t] = sum(ct_impact[t]);  //Cuánto de los ingresos depende de las diferentes variables de control
  }
  
  //-----------------------------------
  //SIGNAL DEFINITION
  //-----------------------------------
  real Y_signal[T];
  
  if(C>0){
    for (t in 1:T){
    Y_signal[t] = Y_temporal[t] + Y_marketing[t] + Y_control[t];  //Ingreso final explicable por el modelo
    }
  }else{
    for (t in 1:T){
    Y_signal[t] = Y_temporal[t] + Y_marketing[t];  //Ingreso final explicable por el modelo  
    }
  }

  //-----------------------------------
  //APPLIED RATIOS FOR CONTROL
  //-----------------------------------
  
  // BIG BLOCKS EFFECT ON ESTIMABLE EFFECTS
  
  real RE_Signal_basalTrend;   // Ratio of Effects over the Y_signal of the beta_0 (intercept sales)
  RE_Signal_basalTrend= ((T*beta_0)+sum(trend))/sum(Y_signal);
    
  real RE_Signal_Y_temporal; // Ratio of Effects over the Y_signal of the Temporal Effects
  RE_Signal_Y_temporal= sum(Y_temporal)/sum(Y_signal);
  
  real RE_Signal_Y_Marketing; // Ratio of Effects over the Y_signal of the Advertising Spends
  RE_Signal_Y_Marketing= sum(Y_marketing)/sum(Y_signal);
  
  real RE_Signal_Y_Control; // Ratio of Effects over the Y_signal of the Advertising Spends
  RE_Signal_Y_Control= sum(Y_control)/sum(Y_signal);
  
  // TEMPORAL EFFECTS
  
  real RE_Temporal_trend; // Ratio of Effects over the Y_temporal of the Trend
  RE_Temporal_trend= sum(trend)/sum(Y_temporal);
  
  real RE_Temporal_seasonal;
  RE_Temporal_seasonal= sum(seasonal)/sum(Y_temporal);
  
  // ADVERTISING EFFECTS
  real RE_advertising_x[M]; // Ratio of Effects over the Y_marketing of the media x[i]
  real RE_Signal_x[M]; // Ratio of Effects over the Y_marketing of the media x[i]
  for (m in 1:M) {
    RE_advertising_x[m] = sum(x_hilladstocked[,m])/sum(Y_marketing);
    RE_Signal_x[m]= sum(x_hilladstocked[,m])/sum(Y_signal);
  }
  
  // CONTROL EFFECTS
  real RE_Signal_ct[C_used]; // Ratio of Effects over the Y_marketing of the media x[i]
    for (c in 1:C_used) {
      RE_Signal_ct[c]= sum(ct_impact[,c])/sum(Y_signal);
      }
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//MODEL
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
model {

    for (t in 1:T){
      if (is_test[t] == 0) {
        Y[t]~normal(Y_signal[t],sigma_e);
      }
    }
  
    //Prioris
  
     beta_0~weibull(beta_0_shape,beta_0_scale); 
     sigma_e~weibull(sigma_e_shape,sigma_e_scale); 
  
     tau~normal(tau_mean,tau_sigma);
     
     psi~normal(0, sigma_psi);
     
     sigma_psi~normal(sigma_psi_mean,sigma_psi_sigma);
  
     
  for (m in 1:M){
       
     theta_x[m]~beta(theta_x_a[m],theta_x_b[m]);
     
     shape_x[m]~weibull(shape_x_shape[m]+1,shape_x_scale[m]);
     scale_x[m]~weibull(scale_x_shape[m],scale_x_scale[m]);
     
     lambda_x[m]~uniform(lambda_x_min[m],lambda_x_max[m]);
     
     beta_x[m]~weibull(beta_x_shape[m],beta_x_scale[m]);   
     gamma_x[m]~weibull(gamma_x_shape[m],gamma_x_scale[m]);
     alpha_x[m]~weibull(alpha_x_shape[m],alpha_x_scale[m]);
  }
  
  for (c in 1:C_used){   
     beta_ct[c]~normal(beta_ct_mu[c],beta_ct_sigma[c]);   
  }

}


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//GENERATED QUANTITIES
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

generated quantities {
  
  // Simulated Synthetic data
  //---------------------------
  real y_sim[T, M_sim];         

  for (t in 1:T) {
    for (m in 1:M_sim) {
      y_sim[t,m] = normal_rng(Y_signal[t], sigma_e);
    }
  }

  // Prediction
  //---------------------------
  vector[T] y_hat;

  for (t in 1:T) {
    y_hat[t] = Y_signal[t];
  }

  // ERROR METRICS on test set
  //----------------------------
    int test_count = 0;
    real<lower=0> y_test_sum = 0;
    real<lower=0> y_test_mean;
    
    if(sum(is_test)>0){
      for (t in 1:T) {
      if (is_test[t] == 1) {
        test_count += 1;
        y_test_sum += Y[t];
      }}
      
      y_test_mean = y_test_sum/test_count;
    }
    
    // MSE - Error Cuadrático Medio ------------------
    real<lower=0> sq_err = 0;
    real<lower=0> MSE;
   
   if(sum(is_test)>0){
    for (t in 1:T) {
      if (is_test[t] == 1) {
        sq_err += square(y_hat[t] - Y[t]);
      }}
      MSE = sq_err / test_count;
    }
    
    // RMSE ------------------
    real<lower=0> RMSE;
    
    if(sum(is_test)>0){
      RMSE = sqrt(MSE);
    }
    // NRMSE ------------------
    real<lower=0> NRMSE;
    
    if(sum(is_test)>0){
      NRMSE= RMSE/y_test_mean;
    }
}
