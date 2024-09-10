//+------------------------------------------------------------------+
//|                                             AllAverages_v4.9.mq4 |
//|                             Copyright © 2007-19, TrendLaboratory |
//+------------------------------------------------------------------+
#property copyright "Copyright © 2007-19, TrendLaboratory"

#property indicator_chart_window
#property indicator_buffers 3
#property indicator_color1  clrTomato
#property indicator_width1  2  
#property indicator_color2  clrDeepSkyBlue
#property indicator_width2  2  
#property indicator_color3  clrDeepSkyBlue
#property indicator_width3  2  

#property strict 


enum ENUM_MA_MODE
{
   SMA,                 // Simple Moving Average
   EMA,                 // Exponential Moving Average
   Wilder,              // Wilder Exponential Moving Average
   LWMA,                // Linear Weighted Moving Average
   SineWMA,             // Sine Weighted Moving Average
   TriMA,               // Triangular Moving Average
   LSMA,                // Least Square Moving Average (or EPMA, Linear Regression Line)
   SMMA,                // Smoothed Moving Average
   HMA,                 // Hull Moving Average by A.Hull
   ZeroLagEMA,          // Zero-Lag Exponential Moving Average
   DEMA,                // Double Exponential Moving Average by P.Mulloy
   T3_basic,            // T3 by T.Tillson (original version)
   ITrend,              // Instantaneous Trendline by J.Ehlers
   Median,              // Moving Median
   GeoMean,             // Geometric Mean
   REMA,                // Regularized EMA by C.Satchwell
   ILRS,                // Integral of Linear Regression Slope
   IE_2,                // Combination of LSMA and ILRS
   TriMAgen,            // Triangular Moving Average generalized by J.Ehlers
   VWMA,                // Volume Weighted Moving Average
   JSmooth,             // M.Jurik's Smoothing
   SMA_eq,              // Simplified SMA
   ALMA,                // Arnaud Legoux Moving Average
   TEMA,                // Triple Exponential Moving Average by P.Mulloy
   T3,                  // T3 by T.Tillson (correct version)
   Laguerre,            // Laguerre filter by J.Ehlers
   MD,                  // McGinley Dynamic
   BF2P,                // Two-pole modified Butterworth filter by J.Ehlers
   BF3P,                // Three-pole modified Butterworth filter by J.Ehlers
   SuperSmu,            // SuperSmoother by J.Ehlers
   Decycler,            // Simple Decycler by J.Ehlers
   eVWMA,               // Modified eVWMA
   EWMA,                // Exponential Weighted Moving Average
   DsEMA,               // Double Smoothed EMA
   TsEMA,               // Triple Smoothed EMA
   VEMA                 // Volume-weighted Exponential Moving Average(V-EMA)
};   


enum ENUM_PRICE
{
   close,               // Close
   open,                // Open
   high,                // High
   low,                 // Low
   median,              // Median
   typical,             // Typical
   weightedClose,       // Weighted Close
   medianBody,          // Median Body (Open+Close)/2
   average,             // Average (High+Low+Open+Close)/4
   trendBiased,         // Trend Biased
   trendBiasedExt,      // Trend Biased(extreme)
   haClose,             // Heiken Ashi Close
   haOpen,              // Heiken Ashi Open
   haHigh,              // Heiken Ashi High   
   haLow,               // Heiken Ashi Low
   haMedian,            // Heiken Ashi Median
   haTypical,           // Heiken Ashi Typical
   haWeighted,          // Heiken Ashi Weighted Close
   haMedianBody,        // Heiken Ashi Median Body
   haAverage,           // Heiken Ashi Average
   haTrendBiased,       // Heiken Ashi Trend Biased
   haTrendBiasedExt     // Heiken Ashi Trend Biased(extreme)   
};

#define pi 3.14159265358979323846

//---- 
input ENUM_TIMEFRAMES   TimeFrame            =           0;       // Timeframe
input ENUM_PRICE        Price                =           0;       // Apply To
input int               MA_Period            =          14;       // Period
input int               MA_Shift             =           0;       // Shift
input ENUM_MA_MODE      MA_Method            =         SMA;       // Method
input bool              ShowInColor          =        true;       // Show In Color
input int               CountBars            =           0;       // Number of bars counted: 0-all bars 

input string            Alerts               = "=== Alerts, Emails & Notifications ===";
input bool              AlertOn              =       false;       // Alert On
input int               AlertShift           =           1;       // Alert Shift:0-current bar,1-previous bar
input int               SoundsNumber         =           5;       // Number of sounds after Signal
input int               SoundsPause          =           5;       // Pause in sec between sounds 
input string            UpTrendSound         = "alert.wav";       // UpTrend Sound
input string            DnTrendSound         = "alert2.wav";      // DownTrend Sound
input bool              EmailOn              =       false;       // Email Alert On
input int               EmailsNumber         =           1;       // Emails Number  
input bool              PushNotificationOn   =       false;       // Push Notification On



double allAvg[];
double uptrend1[];
double uptrend2[];
double iprice[];
double trend[];

//----
int      cBars, timeframe, draw_begin, ma_period, masize; 
string   IndicatorName, TF, short_name, maname;
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
int OnInit()
{
   timeframe = TimeFrame;
   if(timeframe <= Period()) timeframe = Period(); 
   TF = tf(timeframe);
   
   IndicatorDigits(Digits);
//----
   IndicatorBuffers(5); 
   
   SetIndexBuffer(0,  allAvg); SetIndexStyle(0,DRAW_LINE); SetIndexShift(0,MA_Shift*timeframe/Period());
   SetIndexBuffer(1,uptrend1); SetIndexStyle(1,DRAW_LINE); SetIndexShift(1,MA_Shift*timeframe/Period());
   SetIndexBuffer(2,uptrend2); SetIndexStyle(2,DRAW_LINE); SetIndexShift(2,MA_Shift*timeframe/Period());
   SetIndexBuffer(3,  iprice);
   SetIndexBuffer(4,   trend); 
   
   if(MA_Method == EWMA) ma_period = 5*MA_Period; else ma_period = MA_Period;
   if(MA_Method == HMA ) ma_period += (int)MathSqrt(MA_Period);
   if(CountBars == 0) cBars = timeframe/Period()*(iBars(NULL,TimeFrame) - ma_period - 3); else cBars = CountBars*timeframe/Period();
   
   draw_begin = Bars - cBars;
   SetIndexDrawBegin(0,draw_begin);   
   SetIndexDrawBegin(1,draw_begin);
   SetIndexDrawBegin(2,draw_begin);
//---- 
   maname = EnumToString(MA_Method);
   masize = averageSize(MA_Method);
      
   IndicatorName = WindowExpertName(); 
   short_name    = IndicatorName+"["+TF+"]("+maname+","+EnumToString(Price)+","+(string)MA_Period+")";  
   IndicatorShortName(short_name);
   
   SetIndexLabel(0,maname+"("+EnumToString(Price)+","+(string)MA_Period+")");
   SetIndexLabel(1,NULL);
   SetIndexLabel(2,NULL); 
//----   
   if(masize > 0) ArrayResize(tmp,masize);
  
      
   return(INIT_SUCCEEDED);
}
//+------------------------------------------------------------------+
//| AllAverages_v4.9 600+                                            |
//+------------------------------------------------------------------+
int start()
{
   int limit = 0, i, y, shift, counted_bars = IndicatorCounted(); 
   
   if(counted_bars > 0) limit = Bars - counted_bars - 1;
   if(counted_bars < 0) return(0);    
   if(counted_bars < 1)
   {
   limit = MathMin(Bars - 2,(cBars + ma_period)*timeframe/Period());
      for(i=Bars-1;i>=0;i--) 
      {
      allAvg[i]   = EMPTY_VALUE; 
      uptrend1[i] = EMPTY_VALUE; 
      uptrend2[i] = EMPTY_VALUE; 
      }
   }
   
   
//---- 
   if(timeframe != Period())
	{
   limit = MathMax(limit,timeframe/Period());   
      
      for(shift=0;shift<limit;shift++) 
      {	
      y = iBarShift(NULL,TimeFrame,Time[shift]);
            
      allAvg[shift] = iCustom(NULL,TimeFrame,IndicatorName,0,Price,MA_Period,0,MA_Method,ShowInColor,CountBars,
                              "",AlertOn,AlertShift,SoundsNumber,SoundsPause,UpTrendSound,DnTrendSound,EmailOn,EmailsNumber,PushNotificationOn,0,y);    
         
         if(ShowInColor)
         {
         trend[shift] = iCustom(NULL,TimeFrame,IndicatorName,0,Price,MA_Period,0,MA_Method,ShowInColor,CountBars,
                                "",AlertOn,AlertShift,SoundsNumber,SoundsPause,UpTrendSound,DnTrendSound,EmailOn,EmailsNumber,PushNotificationOn,4,y);    
         
      
         uptrend1[shift] = EMPTY_VALUE;
         if(trend[shift] > 0) uptrend1[shift] = allAvg[shift];  
         }
      }
      
      if(ShowInColor)
      {
         for(shift=limit;shift>=0;shift--) 
         {	
         uptrend2[shift] = EMPTY_VALUE;   
            if(allAvg[shift] > allAvg[shift+1]) 
            {
            uptrend2[shift]   = allAvg[shift];
            uptrend2[shift+1] = allAvg[shift+1];
            }	
         }
      }
      
      if(CountBars > 0)
      {
      SetIndexDrawBegin(0,Bars - cBars);   
      SetIndexDrawBegin(1,Bars - cBars);
      SetIndexDrawBegin(2,Bars - cBars);
      }   
   
	return(0);
	}
   
      
   for(shift=limit;shift>=0;shift--) 
   {
      if(Price <= 10) iprice[shift] = getPrice((int)Price,shift);   
      else
      if(Price > 10 && Price <= 21) iprice[shift] = HeikenAshi(0,Price-11,MathMin(Bars-1,cBars+ma_period),shift);
   
   if(shift > MathMin(Bars-ma_period-1,cBars)) continue;
   
   allAvg[shift] = allAveragesOnArray(0,iprice,MA_Period,(int)MA_Method,masize,MathMin(Bars-ma_period-1,cBars),shift);   
              
   trend[shift] = trend[shift+1];
   if(allAvg[shift] > allAvg[shift+1]) trend[shift] = 1;
   if(allAvg[shift] < allAvg[shift+1]) trend[shift] =-1;       
         
      if(ShowInColor)
      {
      uptrend1[shift] = EMPTY_VALUE; 
      uptrend2[shift] = EMPTY_VALUE;
          
         for(i=0;i<2;i++)
         {
            if(trend[shift+i] == 1) 
            {
               if(uptrend1[shift+1+i] == EMPTY_VALUE)
               {
                  if(uptrend1[shift+2+i] == EMPTY_VALUE) 
                  {
                  uptrend1[shift+i]   = allAvg[shift+i];
                  uptrend1[shift+1+i] = allAvg[shift+1+i];
                  uptrend2[shift+i]   = EMPTY_VALUE;
                  }
                  else 
                  {
                  uptrend2[shift+i]   = allAvg[shift+i];
                  uptrend2[shift+1+i] = allAvg[shift+1+i];
                  uptrend1[shift+i]   = EMPTY_VALUE;
                  }
               }
               else
               {
               uptrend1[shift+i] = allAvg[shift+i];
               uptrend2[shift+i] = EMPTY_VALUE;
               }
            }
         }
      }
   }
            
   
   if(AlertOn || EmailOn || PushNotificationOn)
   {
   bool uptrend = trend[AlertShift] > 0 && trend[AlertShift+1] <= 0;                  
   bool dntrend = trend[AlertShift] < 0 && trend[AlertShift+1] >= 0;
         
      if(uptrend || dntrend)
      {
         if(isNewBar(timeframe))
         {
            if(AlertOn)
            {
            BoxAlert(uptrend," : " + maname + "(" + (string)MA_Period + ") is going Up @ "  +DoubleToStr(Close[AlertShift],Digits));   
            BoxAlert(dntrend," : " + maname + "(" + (string)MA_Period + ") is going Down @ "+DoubleToStr(Close[AlertShift],Digits)); 
            }
                   
            if(EmailOn)
            {
            EmailAlert(uptrend,"BUY" ," : " + maname + "(" + (string)MA_Period + ") is going Up @ "  +DoubleToStr(Close[AlertShift],Digits),EmailsNumber); 
            EmailAlert(dntrend,"SELL"," : " + maname + "(" + (string)MA_Period + ") is going Down @ "+DoubleToStr(Close[AlertShift],Digits),EmailsNumber); 
            }
         
            if(PushNotificationOn)
            {
            PushAlert(uptrend," : " + maname + "(" + (string)MA_Period + ") is going Up @ "  +DoubleToStr(Close[AlertShift],Digits));   
            PushAlert(dntrend," : " + maname + "(" + (string)MA_Period + ") is going Down @ "+DoubleToStr(Close[AlertShift],Digits)); 
            }
         }
         else
         {
            if(AlertOn)
            {
            WarningSound(uptrend,SoundsNumber,SoundsPause,UpTrendSound,Time[AlertShift]);
            WarningSound(dntrend,SoundsNumber,SoundsPause,DnTrendSound,Time[AlertShift]);
            }
         }     
      }
   }   
      
      
   if(CountBars > 0)
   {
   SetIndexDrawBegin(0,Bars - cBars);   
   SetIndexDrawBegin(1,Bars - cBars);
   SetIndexDrawBegin(2,Bars - cBars);
   }

//---- 
   return(0);
}

int averageSize(int mode)
{   
   int arraysize;
   
   switch(mode)
   {
   case DEMA     : arraysize = 2; break;
   case T3_basic : arraysize = 6; break;
   case JSmooth  : arraysize = 5; break;
   case TEMA     : arraysize = 4; break;
   case T3       : arraysize = 6; break;
   case Laguerre : arraysize = 4; break;
   case DsEMA    : arraysize = 2; break;
   case TsEMA    : arraysize = 3; break;
   case VEMA     : arraysize = 2; break;
   default       : arraysize = 0; break;
   }
   
   return(arraysize);
}


double   tmp[][1][2], ma[1][4];
datetime prevtime[1];  

double allAveragesOnArray(int index,double& price[],int period,int mode,int arraysize,int cbars,int bar)
{
   int i;
   double MA[4];  
        
   switch(mode) 
	{
	case EMA: case Wilder: case SMMA: case ZeroLagEMA: case DEMA: case T3_basic: case ITrend: case REMA: case JSmooth: case SMA_eq: case TEMA: case T3: case Laguerre: case MD: 
	case BF2P: case BF3P: case SuperSmu: case Decycler: case eVWMA: case DsEMA: case TsEMA: case VEMA:
		
		if(prevtime[index] != Time[bar])
      {
      ma[index][3] = ma[index][2]; 
      ma[index][2] = ma[index][1]; 
      ma[index][1] = ma[index][0]; 
   
      if(arraysize > 0) for(i=0;i<arraysize;i++) tmp[i][index][1] = tmp[i][index][0];
    
      prevtime[index] = Time[bar]; 
      }
   
      if(mode == ITrend || mode == REMA || mode == SMA_eq || (mode >= BF2P && mode < eVWMA)) for(i=0;i<4;i++) MA[i] = ma[index][i]; 
	}
      
   switch(mode)
   {
   case EMA       : ma[index][0] = EMAOnArray(price[bar],ma[index][1],period,cbars,bar); break;
   case Wilder    : ma[index][0] = WilderOnArray(price[bar],ma[index][1],period,cbars,bar); break;  
   case LWMA      : ma[index][0] = LWMAOnArray(price,period,bar); break;
   case SineWMA   : ma[index][0] = SineWMAOnArray(price,period,bar); break;
   case TriMA     : ma[index][0] = TriMAOnArray(price,period,bar); break;
   case LSMA      : ma[index][0] = LSMAOnArray(price,period,bar); break;
   case SMMA      : ma[index][0] = SMMAOnArray(price,ma[index][1],period,cbars,bar); break;
   case HMA       : ma[index][0] = HMAOnArray(price,period,cbars,bar); break;
   case ZeroLagEMA: ma[index][0] = ZeroLagEMAOnArray(price,ma[index][1],period,cbars,bar); break;
   case DEMA      : ma[index][0] = DEMAOnArray(index,0,price[bar],period,1,cbars,bar); break;
   case T3_basic  : ma[index][0] = T3_basicOnArray(index,0,price[bar],period,0.7,cbars,bar); break;
   case ITrend    : ma[index][0] = ITrendOnArray(price,MA,period,cbars,bar); break;
   case Median    : ma[index][0] = MedianOnArray(price,period,cbars,bar); break;
   case GeoMean   : ma[index][0] = GeoMeanOnArray(price,period,cbars,bar); break;
   case REMA      : ma[index][0] = REMAOnArray(price[bar],MA,period,0.5,cbars,bar); break;
   case ILRS      : ma[index][0] = ILRSOnArray(price,period,cbars,bar); break;
   case IE_2      : ma[index][0] = IE2OnArray(price,period,cbars,bar); break;
   case TriMAgen  : ma[index][0] = TriMA_genOnArray(price,period,cbars,bar); break;
   case VWMA      : ma[index][0] = VWMAOnArray(price,period,bar); break;
   case JSmooth   : ma[index][0] = JSmoothOnArray(index,0,price[bar],period,1,cbars,bar); break;
   case SMA_eq    : ma[index][0] = SMA_eqOnArray(price,MA,period,cbars,bar); break;
   case ALMA      : ma[index][0] = ALMAOnArray(price,period,0.85,8,bar); break;
   case TEMA      : ma[index][0] = TEMAOnArray(index,price[bar],period,1,cbars,bar); break;
   case T3        : ma[index][0] = T3OnArray(index,0,price[bar],period,0.7,cbars,bar); break;
   case Laguerre  : ma[index][0] = LaguerreOnArray(index,price[bar],period,4,cbars,bar); break;
   case MD        : ma[index][0] = McGinleyOnArray(price[bar],ma[index][1],period,cbars,bar); break;
   case BF2P      : ma[index][0] = BF2POnArray(price,MA,period,cbars,bar); break;
   case BF3P      : ma[index][0] = BF3POnArray(price,MA,period,cbars,bar); break;
   case SuperSmu  : ma[index][0] = SuperSmuOnArray(price,MA,period,cbars,bar); break;
   case Decycler  : ma[index][0] = DecyclerOnArray(price,MA,period,cbars,bar); return(price[bar] - ma[index][0]); 
   case eVWMA     : ma[index][0] = eVWMAOnArray(price[bar],ma[index][1],period,cbars,bar); break;
   case EWMA      : ma[index][0] = EWMAOnArray(price,period,bar); break;
   case DsEMA     : ma[index][0] = DsEMAOnArray(index,price[bar],period,cbars,bar); break;
   case TsEMA     : ma[index][0] = TsEMAOnArray(index,price[bar],period,cbars,bar); break;
   case VEMA      : ma[index][0] = VEMAOnArray(index,price[bar],period,cbars,bar); break;
   default        : ma[index][0] = SMAOnArray(price,period,bar); break;
   }
   
   return(ma[index][0]);
}

// MA_Method=SMA - Simple Moving Average
double SMAOnArray(double& array[],int per,int bar)
{
   double sum = 0;
   for(int i=0;i<per;i++) sum += array[bar+i];
   
   return(sum/per);
}                
// MA_Method=EMA - Exponential Moving Average
double EMAOnArray(double price,double prev,int per,int cbars,int bar)
{
   double ema;
   
   if(bar >= cbars) ema = price;
   else 
   ema = prev + 2.0/(1 + per)*(price - prev); 
   
   return(ema);
}
// MA_Method=Wilder - Wilder Exponential Moving Average
double WilderOnArray(double price,double prev,int per,int cbars,int bar)
{
   double wilder;
   
   if(bar >= cbars) wilder = price; 
   else 
   wilder = prev + (price - prev)/per; 
   
   return(wilder);
}

// MA_Method=LWMA - Linear Weighted Moving Average 
double LWMAOnArray(double& array[],int per,int bar)
{
   double sum = 0, weight = 0;
   
      for(int i=0;i<per;i++)
      { 
      weight += (per - i);
      sum    += array[bar+i]*(per - i);
      }
   
   if(weight > 0) return(sum/weight); else return(0); 
} 

// MA_Method=SineWMA - Sine Weighted Moving Average
double SineWMAOnArray(double& array[],int per,int bar)
{
   double sum = 0, weight = 0;
  
      for(int i=0;i<per;i++)
      { 
      weight += MathSin(pi*(i + 1)/(per + 1));
      sum    += array[bar+i]*MathSin(pi*(i + 1)/(per + 1)); 
      }
   
   if(weight > 0) return(sum/weight); else return(0); 
}

// MA_Method=TriMA - Triangular Moving Average
double TriMAOnArray(double& array[],int per,int bar)
{
   int len    = (int)MathCeil((per + 1)*0.5);
   double sum = 0;
   
   for(int i=0;i<len;i++) sum += SMAOnArray(array,len,bar+i);
         
   return(sum/len);
}

// MA_Method=LSMA - Least Square Moving Average (or EPMA, Linear Regression Line)
double LSMAOnArray(double& array[],int per,int bar)
{   
   double sum = 0;
   
   for(int i=per;i>=1;i--) sum += (i - (per + 1)/3.0)*array[bar+per-i];
   
   return(sum*6/(per*(per + 1)));
}

// MA_Method=SMMA - Smoothed Moving Average
double SMMAOnArray(double& array[],double prev,int per,int cbars,int bar)
{
   double smma = 0;
   
   if(bar == cbars) smma = SMAOnArray(array,per,bar);
   else 
   if(bar  < cbars)
   {
   double sum = 0;
   for(int i=0;i<per;i++) sum += array[bar+i+1];
   smma = (sum - prev + array[bar])/per;
   }
   
   return(smma);
}                

// MA_Method=HMA - Hull Moving Average by Alan Hull
double HMAOnArray(double& array[],int per,int cbars,int bar)
{
   double hma = 0, _tmp[];
   int len = (int)MathSqrt(per);
   
   ArrayResize(_tmp,len);
   
   if(bar == cbars) hma = array[bar]; 
   else
   if(bar  < cbars)
   {
   for(int i=0;i<len;i++) _tmp[i] = 2*LWMAOnArray(array,per/2,bar+i) - LWMAOnArray(array,per,bar+i);  
   hma = LWMAOnArray(_tmp,len,0); 
   }  

   return(hma);
}

// MA_Method=ZeroLagEMA - Zero-Lag Exponential Moving Average
double ZeroLagEMAOnArray(double& price[],double prev,int per,int cbars,int bar)
{
   int lag = (int)(0.5*(per - 1)); 
   double zema, alpha = 2.0/(1 + (double)per); 
      
   if(bar >= cbars) zema = price[bar];
   else 
   zema = alpha*(2*price[bar] - price[bar+lag]) + (1 - alpha)*prev;
   
   return(zema);
}

// MA_Method=DEMA - Double Exponential Moving Average by Patrick Mulloy
double DEMAOnArray(int index,int num,double price,double per,double v,int cbars,int bar)
{
   double dema = 0, alpha = 2.0/(1 + per);
   
   if(bar == cbars) {dema = price; tmp[num][index][0] = dema; tmp[num+1][index][0] = dema;}
   else 
   if(bar <  cbars) 
   {
   tmp[num  ][index][0] = tmp[num  ][index][1] + alpha*(price              - tmp[num  ][index][1]); 
   tmp[num+1][index][0] = tmp[num+1][index][1] + alpha*(tmp[num][index][0] - tmp[num+1][index][1]); 
   dema                 = tmp[num  ][index][0]*(1+v) - tmp[num+1][index][0]*v;
   }
   
   return(dema);
}

// MA_Method=T3 by T.Tillson
double T3_basicOnArray(int index,int num,double price,int per,double v,int cbars,int bar)
{
   double dema1, dema2, T3 = 0;
   
   if(bar == cbars) 
   {
   T3 = price; 
   for(int k=0;k<6;k++) tmp[num+k][index][0] = price;
   }
   else 
   if(bar < cbars) 
   {
   dema1 = DEMAOnArray(index,num  ,price,per,v,cbars,bar); 
   dema2 = DEMAOnArray(index,num+2,dema1,per,v,cbars,bar); 
   T3    = DEMAOnArray(index,num+4,dema2,per,v,cbars,bar);
   }
   
   return(T3);
}

// MA_Method=ITrend - Instantaneous Trendline by J.Ehlers
double ITrendOnArray(double& price[],double& array[],int per,int cbars,int bar)
{
   double it = 0, alpha = 2.0/(per + 1);
   
   if(bar < cbars && array[1] > 0 && array[2] > 0) it = (alpha - 0.25*alpha*alpha)*price[bar] + 0.5*alpha*alpha*price[bar+1] 
                                                       -(alpha - 0.75*alpha*alpha)*price[bar+2] + 2*(1 - alpha)*array[1] 
                                                       -(1 - alpha)*(1 - alpha)*array[2];
   else it = (price[bar] + 2*price[bar+1] + price[bar+2])/4;
      
   return(it);
}

// MA_Method=Median - Moving Median
double MedianOnArray(double& price[],int per,int cbars,int bar)
{
   double median = 0, array[];
   ArrayResize(array,per);
   
   if(bar <= cbars)
   {
   for(int i=0;i<per;i++) array[i] = price[bar+i];
   ArraySort(array,WHOLE_ARRAY,0,MODE_DESCEND);
   
   int num = (int)MathRound((per - 1)*0.5); 
   if(MathMod(per,2) > 0) median = array[num]; else median = 0.5*(array[num] + array[num+1]);
   }
    
   return(median); 
}

// MA_Method=GeoMean - Geometric Mean
double GeoMeanOnArray(double& price[],int per,int cbars,int bar)
{
   double gmean = 0;
   
   if(bar < cbars)
   { 
   gmean = MathPow(price[bar],1.0/per); 
   for(int i=1;i<per;i++) gmean *= MathPow(price[bar+i],1.0/per); 
   }
   else 
   if(bar == cbars) gmean = SMAOnArray(price,per,bar);
   
   return(gmean);
}

// MA_Method=REMA - Regularized EMA by Chris Satchwell 
double REMAOnArray(double price,double& array[],int per,double lambda,int cbars,int bar)
{
   double rema, alpha =  2.0/(per + 1);
   
   if(bar < cbars  && array[1] > 0 && array[2] > 0) rema = (array[1]*(1 + 2*lambda) + alpha*(price - array[1]) - lambda*array[2])/(1 + lambda); 
   else rema = price; 
   
   return(rema);
}
// MA_Method=ILRS - Integral of Linear Regression Slope 
double ILRSOnArray(double& price[],int per,int cbars,int bar)
{
   double slope, sum  = per*(per - 1)*0.5;
   double ilrs = 0, sum2 = (per - 1)*per*(2*per - 1)/6.0;
   
   if(bar < cbars)
   {  
   double sum1 = 0;
   double sumy = 0;
      for(int i=0;i<per;i++)
      { 
      sum1 += i*price[bar+i];
      sumy += price[bar+i];
      }
   double num1 = per*sum1 - sum*sumy;
   double num2 = sum*sum - per*sum2;
   
   if(num2 != 0) slope = num1/num2; else slope = 0; 
   ilrs = slope + SMAOnArray(price,per,bar);
   }
   
   return(ilrs);
}
// MA_Method=IE/2 - Combination of LSMA and ILRS 
double IE2OnArray(double& price[],int per,int cbars,int bar)
{
   double ie = 0;
   
   if(bar < cbars) ie = 0.5*(ILRSOnArray(price,per,cbars,bar) + LSMAOnArray(price,per,bar));
      
   return(ie); 
}
 
// MA_Method=TriMAgen - Triangular Moving Average Generalized by J.Ehlers
double TriMA_genOnArray(double& array[],int per,int cbars,int bar)
{
   int len1 = (int)MathFloor((per + 1)*0.5);
   int len2 = (int)MathCeil ((per + 1)*0.5);
   double sum = 0;
   
   if(bar < cbars)
       for(int i = 0;i < len2;i++) sum += SMAOnArray(array,len1,bar+i);
   
   return(sum/len2);
}

// MA_Method=VWMA - Volume Weighted Moving Average 
double VWMAOnArray(double& array[],int per,int bar)
{
   double sum = 0, weight = 0, maxvol = 0, vwma = 0;
   
      for(int i=0;i<per;i++)
      { 
      weight += (double)Volume[bar+i];
      sum    += array[bar+i]*Volume[bar+i];
      }
     
   if(weight > 0) return(sum/weight); else return(0); 
} 

// MA_Method=JSmooth - Smoothing by Mark Jurik
double JSmoothOnArray(int index,int num,double price,int per,double power,int cbars,int bar)
{
   double beta  = 0.45*(per - 1)/(0.45*(per - 1) + 2);
	double alpha = MathPow(beta,power);
	
	if(bar == cbars) {tmp[num+4][index][0] = price; tmp[num+0][index][0] = price; tmp[num+2][index][0] = price;}
	else 
   if(bar <  cbars) 
   {
	tmp[num+0][index][0] = (1 - alpha)*price + alpha*tmp[num+0][index][1];
	tmp[num+1][index][0] = (price - tmp[num+0][index][0])*(1-beta) + beta*tmp[num+1][index][1];
	tmp[num+2][index][0] = tmp[num+0][index][0] + tmp[num+1][index][0];
	tmp[num+3][index][0] = (tmp[num+2][index][0] - tmp[num+4][index][1])*MathPow((1-alpha),2) + MathPow(alpha,2)*tmp[num+3][index][1];
	tmp[num+4][index][0] = tmp[num+4][index][1] + tmp[num+3][index][0]; 
   }
   return(tmp[num+4][index][0]);
}

// MA_Method=SMA_eq     - Simplified SMA
double SMA_eqOnArray(double& price[],double& array[],int per,int cbars,int bar)
{
   double sma = 0;
   
   if(bar == cbars) sma = SMAOnArray(price,per,bar);
   else 
   if(bar <  cbars) sma = (price[bar] - price[bar+per])/per + array[1]; 
   
   return(sma);
}                        		

// MA_Method=ALMA by Arnaud Legoux / Dimitris Kouzis-Loukas / Anthony Cascino
double ALMAOnArray(double& price[],int per,double offset,double sigma,int bar)
{
   double m = MathFloor(offset*(per - 1)), s = per/sigma, w, sum = 0, wsum = 0;		
	
	for(int i=0;i<per;i++) 
	{
	w     = MathExp(-((i - m)*(i - m))/(2*s*s));
   wsum += w;
   sum  += price[bar+(per-1-i)]*w; 
   }
   
   if(wsum != 0) return(sum/wsum); else return(0);
}   

// MA_Method=TEMA - Triple Exponential Moving Average by Patrick Mulloy
double TEMAOnArray(int index,double price,int per,double v,int cbars,int bar)
{
   double alpha = 2.0/(per+1);
	
	if(bar == cbars) {tmp[0][index][0] = price; tmp[1][index][0] = price; tmp[2][index][0] = price;}
	else 
   if(bar <  cbars) 
   {
	tmp[0][index][0] = tmp[0][index][1] + alpha *(price            - tmp[0][index][1]);
	tmp[1][index][0] = tmp[1][index][1] + alpha *(tmp[0][index][0] - tmp[1][index][1]);
	tmp[2][index][0] = tmp[2][index][1] + alpha *(tmp[1][index][0] - tmp[2][index][1]);
	tmp[3][index][0] = tmp[0][index][0] + v*(tmp[0][index][0] + v*(tmp[0][index][0]-tmp[1][index][0]) - tmp[1][index][0] - v*(tmp[1][index][0] - tmp[2][index][0])); 
	}
   
   return(tmp[3][index][0]);
}

// MA_Method=T3 by T.Tillson (correct version) 
double T3OnArray(int index,int num,double price,int per,double v,int cbars,int bar)
{
   double len = MathMax((per + 5.0)/3.0 - 1,1), dema1, dema2, T3 = 0;
   
   if(bar == cbars ) 
   {
   T3 = price; 
   for(int k=0;k<6;k++) tmp[num+k][index][0] = T3;
   }
   else 
   if(bar < cbars) 
   {
   dema1 = DEMAOnArray(index,num  ,price,len,v,cbars,bar); 
   dema2 = DEMAOnArray(index,num+2,dema1,len,v,cbars,bar); 
   T3    = DEMAOnArray(index,num+4,dema2,len,v,cbars,bar);
   }
      
   return(T3);
}

// MA_Method=Laguerre filter by J.Ehlers
double LaguerreOnArray(int index,double price,int per,int order,int cbars,int bar)
{
   double gamma = 1 - 10.0/(per + 9);
   double aPrice[];
   
   ArrayResize(aPrice,order);
   
   for(int i=0;i<order;i++)
   {
      if(bar >= cbars) tmp[i][index][0] = price;
      else
      {
         if(i == 0) tmp[i][index][0] = (1 - gamma)*price + gamma*tmp[i][index][1];
         else
         tmp[i][index][0] = -gamma * tmp[i-1][index][0] + tmp[i-1][index][1] + gamma * tmp[i][index][1];
      
      aPrice[i] = tmp[i][index][0];
      }
   }
   double laguerre = TriMA_genOnArray(aPrice,order,cbars,0);  

   return(laguerre);
}

// MA_Method=MD - McGinley Dynamic
double McGinleyOnArray(double price,double prev,int per,int cbars,int bar)
{
   double md = 0;
   
   if(bar == cbars) md = price;
   else 
   if(bar <  cbars) 
      if(prev != 0) md = prev + (price - prev)/(per*MathPow(price/prev,4)/2); 
      else md = price;

   return(md);
}

// MA_Method=BF2P - Two-pole modified Butterworth filter
double BF2POnArray(double& price[],double& array[],int per,int cbars,int bar)
{
   double a  = MathExp(-1.414*pi/per);
   double b  = 2*a*MathCos(1.414*1.25*pi/per);
   double c2 = b;
   double c3 = -a*a;
   double c1 = 1 - c2 - c3, bf2p;
   
   if(bar < cbars && array[1] > 0 && array[2] > 0) {bf2p = c1*(price[bar] + 2*price[bar+1] + price[bar+2])/4 + c2*array[1] + c3*array[2];
   if(bar == cbars-1) Print("bar=",bar," cbars=",cbars," bf2p=",bf2p," array[1]=",array[1]," array[2]=",array[2]);
   }
   else bf2p = (price[bar] + 2*price[bar+1] + price[bar+2])/4;
   
   return(bf2p);
}

// MA_Method=BF3P - Three-pole modified Butterworth filter
double BF3POnArray(double& price[],double& array[],int per,int cbars,int bar)
{
   double a  = MathExp(-pi/per);
   double b  = 2*a*MathCos(1.738*pi/per);
   double c  = a*a;
   double d2 = b + c;
   double d3 = -(c + b*c);
   double d4 = c*c;
   double d1 = 1 - d2 - d3 - d4, bf3p;
   
   if(bar < cbars && array[1] > 0 && array[2] > 0 && array[3] > 0) bf3p = d1*(price[bar] + 3*price[bar+1] + 3*price[bar+2] + price[bar+3])/8 + d2*array[1] + d3*array[2] + d4*array[3];
   else bf3p = (price[bar] + 3*price[bar+1] + 3*price[bar+2] + price[bar+3])/8;
   
   return(bf3p);
}

// MA_Method=SuperSmu - SuperSmoother filter
double SuperSmuOnArray(double& price[],double& array[],int per,int cbars,int bar)
{
   double a  = MathExp(-1.414*pi/per);
   double b  = 2*a*MathCos(1.414*pi/per);
   double c2 = b;
   double c3 = -a*a;
   double c1 = 1 - c2 - c3, supsm;
   
   if(bar < cbars && array[1] > 0 && array[2] > 0) supsm = c1*(price[bar] + price[bar+1])/2 + c2*array[1] + c3*array[2]; else supsm = (price[bar] + price[bar+1])/2;
   
   return(supsm);
}

// MA_Method=Decycler - Simple Decycler by J.Ehlers
double DecyclerOnArray(double& price[],double& hp[],int per,int cbars,int bar)
{
   double alpha1 = (MathCos(1.414*pi/per) + MathSin(1.414*pi/per) - 1)/MathCos(1.414*pi/per);
  
   if(bar > cbars - 4) return(0);
   
   hp[0] = (1 - alpha1/2)*(1 - alpha1/2)*(price[bar] - 2*price[bar+1] + price[bar+2]) + 2*(1 - alpha1)*hp[1] - (1 - alpha1)*(1 - alpha1)*hp[2];		
       
   return(hp[0]);
}

// MA_Method=eVWMA - Elastic Volume Weighted Moving Average by C.Fries
double eVWMAOnArray(double price,double prev,int per,int cbars,int bar)
{
   double evwma;
 
   if(bar < cbars && prev > 0)
   {
   double max = 0;
   for(int i=0;i<per;i++) max = MathMax(max,Volume[bar+i]);
      
   double diff = 3*max - Volume[bar];
          
      if(diff < 0) evwma = prev;
      else 
      evwma = (diff*prev + Volume[bar]*price)/(3*max); 
   }
   else evwma = price;
    
   return(evwma);
}

// MA_Method=EWMA - Exponential Weighted Moving Average 
double EWMAOnArray(double& array[],int per,int bar)
{
   double sum = 0, weight = 0, alpha = 2.0/(1 + per);
   
      for(int i=0;i<ma_period;i++)
      { 
      weight += alpha*MathPow(1-alpha,i);
      sum    += array[bar+i]*alpha*MathPow(1-alpha,i);
      }
   
   if(weight > 0) return(sum/weight); else return(0); 
} 

// MA_Method=DsEMA - Double Smoothed EMA
double DsEMAOnArray(int index,double price,int per,int cbars,int bar)
{
   double alpha = 4.0/(per + 3);
   
   if(bar == cbars) 
   {
   for(int k=0;k<2;k++) tmp[k][index][0] = price;
   return(price);
   }
   else 
   if(bar < cbars) 
   {
   tmp[0][index][0] = tmp[0][index][1] + alpha*(price            - tmp[0][index][1]);
	tmp[1][index][0] = tmp[1][index][1] + alpha*(tmp[0][index][0] - tmp[1][index][1]);
   }
      
   return(tmp[1][index][0]);
}

// MA_Method=TsEMA - Triple Smoothed EMA
double TsEMAOnArray(int index,double price,int per,int cbars,int bar)
{
   double alpha = 6.0/(per + 5);
   
   if(bar == cbars) 
   {
   for(int k=0;k<3;k++) tmp[k][index][0] = price;
   return(price);
   }
   else 
   if(bar < cbars) 
   {
   tmp[0][index][0] = tmp[0][index][1] + alpha*(price            - tmp[0][index][1]);
	tmp[1][index][0] = tmp[1][index][1] + alpha*(tmp[0][index][0] - tmp[1][index][1]);
   tmp[2][index][0] = tmp[2][index][1] + alpha*(tmp[1][index][0] - tmp[2][index][1]);
   }
      
   return(tmp[2][index][0]);
}

// MA_Method=VEMA - Volume-weighted Exponential Moving Average(V-EMA)
double VEMAOnArray(int index,double price,int per,int cbars,int bar)
{
   double vema = price, alpha = 2.0/(per + 1);
   
   if(bar == cbars) 
   {
   tmp[0][index][0] = (1 - alpha)*price*Volume[bar];
   tmp[1][index][0] = (1 - alpha)*Volume[bar];
   }
   else 
   if(bar < cbars) 
   {
   tmp[0][index][0] = tmp[0][index][1] + alpha*(price*Volume[bar] - tmp[0][index][1]);
	tmp[1][index][0] = tmp[1][index][1] + alpha*(Volume[bar]       - tmp[1][index][1]);
   }
  
   if(tmp[1][index][0] > 0) vema = tmp[0][index][0]/tmp[1][index][0];
   
   return(vema);
}


double getPrice(int price,int bar)
{
   double close = Close[bar];
   double open  = Open [bar];
   double high  = High [bar];
   double low   = Low  [bar];
   
   switch(price)
   {   
   case  0: return(close); break;
   case  1: return(open ); break;
   case  2: return(high ); break;
   case  3: return(low  ); break;
   case  4: return((high + low)/2); break;
   case  5: return((high + low + close)/3); break;
   case  6: return((high + low + 2*close)/4); break;
   case  7: return((close + open)/2); break;
   case  8: return((high + low + close + open)/4); break;
   case  9: if(close > open) return((high + close)/2); else return((low + close)/2); break;
   case 10: if(close > open) return(high); else return(low); break;
   default: return(close); break;
   }
}

// HeikenAshi Price
double   haClose[1][2], haOpen[1][2], haHigh[1][2], haLow[1][2];
datetime prevhatime[1];

double HeikenAshi(int index,int price,int cbars,int bar)
{ 
   if(prevhatime[index] != Time[bar])
   {
   haClose[index][1] = haClose[index][0];
   haOpen [index][1] = haOpen [index][0];
   haHigh [index][1] = haHigh [index][0];
   haLow  [index][1] = haLow  [index][0];
   prevhatime[index] = Time[bar];
   }
   
   if(bar == cbars) 
   {
   haClose[index][0] = Close[bar];
   haOpen [index][0] = Open [bar];
   haHigh [index][0] = High [bar];
   haLow  [index][0] = Low  [bar];
   }
   else
   if(bar <  cbars)
   {
   haClose[index][0] = (Open[bar] + High[bar] + Low[bar] + Close[bar])/4;
   haOpen [index][0] = (haOpen[index][1] + haClose[index][1])/2;
   haHigh [index][0] = MathMax(High[bar],MathMax(haOpen[index][0],haClose[index][0]));
   haLow  [index][0] = MathMin(Low [bar],MathMin(haOpen[index][0],haClose[index][0]));
   }
   
   switch(price)
   {
   case  0: return(haClose[index][0]); break;
   case  1: return(haOpen [index][0]); break;
   case  2: return(haHigh [index][0]); break;
   case  3: return(haLow  [index][0]); break;
   case  4: return((haHigh [index][0] + haLow [index][0])/2); break;
   case  5: return((haHigh[index][0] + haLow[index][0] +   haClose[index][0])/3); break;
   case  6: return((haHigh[index][0] + haLow[index][0] + 2*haClose[index][0])/4); break;
   case  7: return((haClose[index][0] + haOpen[index][0])/2); break;
   case  8: return((haHigh[index][0] + haLow[index][0] +   haClose[index][0] + haOpen[index][0])/4); break;
   case  9: if(haClose[index][0] > haOpen[index][0]) return((haHigh[index][0] + haClose[index][0])/2); else return((haLow[index][0] + haClose[index][0])/2); break;
   case 10: if(haClose[index][0] > haOpen[index][0]) return(haHigh[index][0]); else return(haLow[index][0]); break;
   default: return(haClose[index][0]); break;
   }
} 

datetime prevnbtime;

bool isNewBar(int tf)
{
   bool res = false;
   
   if(tf >= 0)
   {
      if(iTime(NULL,tf,0) != prevnbtime)
      {
      res   = true;
      prevnbtime = iTime(NULL,tf,0);
      }   
   }
   else res = true;
   
   return(res);
}

string prevmess;
 
bool BoxAlert(bool cond,string text)   
{      
   string mess = IndicatorName + "("+Symbol()+","+TF + ")" + text;
   
   if (cond && mess != prevmess)
	{
	Alert (mess);
	prevmess = mess; 
	return(true);
	} 
  
   return(false);  
}

datetime pausetime;

bool Pause(int sec)
{
   if(TimeCurrent() >= pausetime + sec) {pausetime = TimeCurrent(); return(true);}
   
   return(false);
}

datetime warningtime;

void WarningSound(bool cond,int num,int sec,string sound,datetime curtime)
{
   static int i;
   
   if(cond)
   {
   if(curtime != warningtime) i = 0; 
   if(i < num && Pause(sec)) {PlaySound(sound); warningtime = curtime; i++;}       	
   }
}

string prevemail;

bool EmailAlert(bool cond,string text1,string text2,int num)   
{      
   string subj = "New " + text1 +" Signal from " + IndicatorName + "!!!";    
   string mess = IndicatorName + "("+Symbol()+","+TF + ")" + text2;
   
   if (cond && mess != prevemail)
	{
	if(subj != "" && mess != "") for(int i=0;i<num;i++) SendMail(subj, mess);  
	prevemail = mess; 
	return(true);
	} 
  
   return(false);  
}

string prevpush;
 
bool PushAlert(bool cond,string text)   
{      
   string push = IndicatorName + "("+Symbol() + "," + TF + ")" + text;
   
   if(cond && push != prevpush)
	{
	SendNotification(push);
	
	prevpush = push; 
	return(true);
	} 
  
   return(false);  
}


string tf(int itimeframe)
{
   string result = "";
   
   switch(itimeframe)
   {
   case PERIOD_M1:   result = "M1" ;
   case PERIOD_M5:   result = "M5" ;
   case PERIOD_M15:  result = "M15";
   case PERIOD_M30:  result = "M30";
   case PERIOD_H1:   result = "H1" ;
   case PERIOD_H4:   result = "H4" ;
   case PERIOD_D1:   result = "D1" ;
   case PERIOD_W1:   result = "W1" ;
   case PERIOD_MN1:  result = "MN1";
   default:          result = "N/A";
   }
   
   if(result == "N/A")
   {
   if(itimeframe <  PERIOD_H1 ) result = "M"  + (string)itimeframe;
   if(itimeframe >= PERIOD_H1 ) result = "H"  + (string)(itimeframe/PERIOD_H1);
   if(itimeframe >= PERIOD_D1 ) result = "D"  + (string)(itimeframe/PERIOD_D1);
   if(itimeframe >= PERIOD_W1 ) result = "W"  + (string)(itimeframe/PERIOD_W1);
   if(itimeframe >= PERIOD_MN1) result = "MN" + (string)(itimeframe/PERIOD_MN1);
   }
   
   return(result); 
}
