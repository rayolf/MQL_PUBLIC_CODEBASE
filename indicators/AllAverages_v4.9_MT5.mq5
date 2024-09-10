//+------------------------------------------------------------------+
//|                                             AllAverages_v4.9.mq5 |
//|                                Copyright © 2019, TrendLaboratory |
//|            http://finance.groups.yahoo.com/group/TrendLaboratory |
//|                                   E-mail: igorad2003@yahoo.co.uk |
//+------------------------------------------------------------------+

#property copyright "Copyright © 2019, TrendLaboratory"
#property link      "http://finance.groups.yahoo.com/group/TrendLaboratory"
#property link      "http://newdigital-world.com/forum.php"



#property indicator_chart_window
#property indicator_buffers 3
#property indicator_plots   1

#property indicator_type1   DRAW_COLOR_LINE
#property indicator_color1  clrYellow,clrDeepSkyBlue,clrOrangeRed
#property indicator_width1  2

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
   HMA,                 // Hull Moving Average by Alan Hull
   ZeroLagEMA,          // Zero-Lag Exponential Moving Average
   DEMA,                // Double Exponential Moving Average by Patrick Mulloy
   T3_basic,            // T3 by T.Tillson (original version)
   ITrend,              // Instantaneous Trendline by J.Ehlers
   Median,              // Moving Median
   GeoMean,             // Geometric Mean
   REMA,                // Regularized EMA by Chris Satchwell
   ILRS,                // Integral of Linear Regression Slope
   IE_2,                // Combination of LSMA and ILRS
   TriMAgen,            // Triangular Moving Average generalized by J.Ehlers
   VWMA,                // Volume Weighted Moving Average
   JSmooth,             // Smoothing by Mark Jurik
   SMA_eq,              // Simplified SMA
   ALMA,                // Arnaud Legoux Moving Average
   TEMA,                // Triple Exponential Moving Average by Patrick Mulloy
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


double  moving[];
double  trend[];
double  iprice[];


ENUM_TIMEFRAMES  tf;

int      mtf_handle, draw_begin, masize, cBars, ma_period;
string   TF, IndicatorName, short_name, maname;
double   mtf_moving[1], mtf_trend[1];
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
int OnInit()
{
   if(TimeFrame <= _Period) tf = _Period; else tf = TimeFrame; 
   TF = timeframeToString(tf);   
   
   IndicatorSetInteger(INDICATOR_DIGITS,_Digits);

//--- indicator buffers mapping 
   SetIndexBuffer(0,moving,        INDICATOR_DATA); PlotIndexSetInteger(0,PLOT_COLOR_INDEXES,3);
   SetIndexBuffer(1, trend, INDICATOR_COLOR_INDEX);
   SetIndexBuffer(2,iprice,INDICATOR_CALCULATIONS);
   
//---   
   PlotIndexSetInteger(0,PLOT_SHIFT,MA_Shift*PeriodSeconds(tf)/PeriodSeconds(Period())); 
   PlotIndexSetInteger(1,PLOT_SHIFT,MA_Shift*PeriodSeconds(tf)/PeriodSeconds(Period())); 
//--- 
   maname = EnumToString(MA_Method);
   masize = averageSize(MA_Method);
   
   IndicatorName = MQLInfoString((int)MQL5_PROGRAM_NAME);
   short_name    = IndicatorName+"["+TF+"]("+EnumToString(Price)+","+(string)MA_Period+","+maname+")";
   IndicatorSetString(INDICATOR_SHORTNAME,short_name);
//---
   if(MA_Method == EWMA) ma_period = 5*MA_Period; else ma_period = MA_Period;
   if(MA_Method == HMA ) ma_period += (int)MathSqrt(MA_Period);
   
   ArrayResize(tmp,masize);
   
   if(TimeFrame > 0 && tf > _Period) 
   mtf_handle = iCustom(Symbol(),tf,IndicatorName,tf,Price,MA_Period,0,MA_Method,ShowInColor,CountBars,      
                        "",AlertOn,AlertShift,SoundsNumber,SoundsPause,UpTrendSound,DnTrendSound,EmailOn,EmailsNumber,PushNotificationOn);

//--- initialization done
   return(INIT_SUCCEEDED);
}
//+------------------------------------------------------------------+
//| Custom indicator iteration function                              |
//+------------------------------------------------------------------+
int OnCalculate(const int rates_total,const int prev_calculated,
                const datetime &Time[],
                const double   &Open[],
                const double   &High[],
                const double   &Low[],
                const double   &Close[],
                const long     &TickVolume[],
                const long     &Volume[],
                const int      &Spread[])
{
   int x, y, shift, limit, mtflimit, copied = 0;
   datetime mtf_time;
   
//--- preliminary calculations
   if(prev_calculated == 0) 
   {
   if(CountBars == 0) cBars = rates_total - ma_period; else cBars = CountBars*PeriodSeconds(tf)/PeriodSeconds(Period());
   limit    = rates_total - cBars - ma_period; 
   mtflimit = rates_total - 1;
   ArrayInitialize(moving,EMPTY_VALUE);
   ArrayInitialize(trend ,0);
   }
   else 
   {
   limit    = rates_total - 1;
   mtflimit = PeriodSeconds(tf)/PeriodSeconds(Period()) + 1;
   }    
 
   
//--- the main loop of calculations
   if(tf > Period())
   { 
   ArraySetAsSeries(Time,true);   
  
      for(shift=0,y=0;shift<mtflimit;shift++)
      {
      if(Time[shift] < iTime(NULL,TimeFrame,y)) y++; 
      mtf_time = iTime(NULL,TimeFrame,y);
      
      copied = CopyBuffer(mtf_handle,0,mtf_time,mtf_time,mtf_moving);
      if(copied <= 0) return(0);
      x = rates_total - shift - 1;
      moving[x] = mtf_moving[0];
           
         if(ShowInColor)
         {
         copied = CopyBuffer(mtf_handle,1,mtf_time,mtf_time,mtf_trend);   
         if(copied <= 0) return(0);
         trend[x] = mtf_trend[0];   
         }
         else trend[x] = 0; 
      }
   
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,rates_total - cBars);   
   }
   else
   for(shift=limit;shift<rates_total;shift++)
   {
      if(Price <= 10) iprice[shift] = getPrice(Price,Close[shift],Open[shift],High[shift],Low[shift],shift);  
      else
      if(Price > 10 && Price <= 21) iprice[shift] = HeikenAshi(0,(int)Price-11,Close[shift],Open[shift],High[shift],Low[shift],rates_total-cBars-ma_period,shift);
   
           
   int startbar = rates_total - cBars;
  
   if(shift < startbar - ma_period + 1) continue;
      
   moving[shift] = allAveragesOnArray(0,iprice,MA_Period,MA_Method,masize,startbar,rates_total,shift);  
    
      if(shift > 0)
      {
         if(ShowInColor && moving[shift-1] != EMPTY_VALUE)
         {
         trend[shift] = trend[shift-1];
         if(moving[shift] > moving[shift-1]) trend[shift] = 1;
         if(moving[shift] < moving[shift-1]) trend[shift] = 2;   
         }    
         else trend[shift] = 0; 
      }
   }
   
   
   if(AlertOn || EmailOn || PushNotificationOn)
   {
   int alertbar = rates_total - 1 - AlertShift;
   bool uptrend = trend[alertbar] == 1 && trend[alertbar-1] == 2;                  
   bool dntrend = trend[alertbar] == 2 && trend[alertbar-1] == 1;
   
      if(uptrend || dntrend)
      {
         if(isNewBar(tf))
         {
            if(AlertOn)
            {
            BoxAlert(uptrend," : " + maname + "(" + (string)MA_Period + ") is going Up @ "  +DoubleToString(Close[alertbar],_Digits));   
            BoxAlert(dntrend," : " + maname + "(" + (string)MA_Period + ") is going Down @ "+DoubleToString(Close[alertbar],_Digits)); 
            }
                   
            if(EmailOn)
            {
            EmailAlert(uptrend,"BUY" ," : " + maname + "(" + (string)MA_Period + ") is going Up @ "  +DoubleToString(Close[alertbar],_Digits),EmailsNumber); 
            EmailAlert(dntrend,"SELL"," : " + maname + "(" + (string)MA_Period + ") is going Down @ "+DoubleToString(Close[alertbar],_Digits),EmailsNumber); 
            }
         
            if(PushNotificationOn)
            {
            PushAlert(uptrend," : " + maname + "(" + (string)MA_Period + ") is going Up @ "  +DoubleToString(Close[alertbar],_Digits));   
            PushAlert(dntrend," : " + maname + "(" + (string)MA_Period + ") is going Down @ "+DoubleToString(Close[alertbar],_Digits)); 
            }
         }
         else
         {
            if(AlertOn)
            {
            WarningSound(uptrend,SoundsNumber,SoundsPause,UpTrendSound,Time[alertbar]);
            WarningSound(dntrend,SoundsNumber,SoundsPause,DnTrendSound,Time[alertbar]);
            }
         }   
      }
   }
   
   
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,rates_total - cBars);  
//--- done       
   return(rates_total);
}
//+------------------------------------------------------------------+
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
int      prevbar[1];  

double allAveragesOnArray(int index,const double& price[],int period,int mode,int arraysize,int cbars,int rates,int bar)
{
   int i;
   
   if(period == 1) return(price[bar]);
      
   double MA[4];  
        
   switch(mode) 
	{
	case EMA: case Wilder: case SMMA: case ZeroLagEMA: case DEMA: case T3_basic: case ITrend: case REMA: case JSmooth: case SMA_eq: case TEMA: case T3: case Laguerre: case MD: 
	case BF2P: case BF3P: case SuperSmu: case Decycler: case eVWMA: case DsEMA: case TsEMA: case VEMA:
	
      if(prevbar[index] != bar)
      {
      ma[index][3] = ma[index][2]; 
      ma[index][2] = ma[index][1]; 
      ma[index][1] = ma[index][0]; 
   
      if(arraysize > 0) for(i=0;i<arraysize;i++) tmp[i][index][1] = tmp[i][index][0];
       
      prevbar[index] = bar; 
      }
   
      if(mode == ITrend || mode == REMA || mode == SMA_eq || (mode >= BF2P && mode < eVWMA)) for(i=0;i<4;i++) MA[i] = ma[index][i]; 
   }
   
   switch(mode)
   {
   case EMA       : ma[index][0] = EMAOnArray(price[bar],ma[index][1],period,cbars,bar); break;
   case Wilder    : ma[index][0] = WilderOnArray(price[bar],ma[index][1],period,cbars,bar); break;  
   case LWMA      : ma[index][0] = LWMAOnArray(0,price,period,bar); break;
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
   case TriMAgen  : ma[index][0] = TriMA_genOnArray(0,price,period,cbars,bar); break;
   case VWMA      : ma[index][0] = VWMAOnArray(price,period,rates,bar); break;
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
   case eVWMA     : ma[index][0] = eVWMAOnArray(price[bar],ma[index][1],period,cbars,rates,bar); break;
   case EWMA      : ma[index][0] = EWMAOnArray(price,period,bar); break;
   case DsEMA     : ma[index][0] = DsEMAOnArray(index,price[bar],period,cbars,bar); break;
   case TsEMA     : ma[index][0] = TsEMAOnArray(index,price[bar],period,cbars,bar); break;
   case VEMA      : ma[index][0] = VEMAOnArray(index,price[bar],period,cbars,rates,bar); break;
   default        : ma[index][0] = SMAOnArray(0,price,period,bar); break;
   }
   
   return(ma[index][0]);
}

// MA_Method=0: SMA - Simple Moving Average
double SMAOnArray(int mode,const double& array[],int per,int bar)
{
   double sum = 0;
   
   if(mode == 0)
   {
      if(bar >= per) for(int i=0;i<per;i++) sum += array[bar-i];
   }
   else for(int i=0;i<per;i++) sum += array[bar+i];
   
   return(sum/per);
}
                           
// MA_Method=1: EMA - Exponential Moving Average
double EMAOnArray(const double price,double prev,int per,int cbars,int bar)
{
   double ema;
   
   if(bar < cbars) ema = price; else ema = prev + 2.0/(1 + per)*(price - prev); 
   
   return(ema);
}

// MA_Method=2: Wilder - Wilder Exponential Moving Average
double WilderOnArray(double price,double prev,int per,int cbars,int bar)
{
   double wilder;
   
   if(bar < cbars) wilder = price; else wilder = prev + (price - prev)/per; 
   
   return(wilder);
}

// MA_Method=3: LWMA - Linear Weighted Moving Average 
double LWMAOnArray(int mode,const double& array[],int per,int bar)
{
   double lwma = EMPTY_VALUE, sum = 0, weight = 0;
   
   if(mode == 0)
   {
      if(bar > per)
      {
         for(int i=0;i<per;i++)
         { 
         weight += (per - i);
         sum    += array[bar-i]*(per - i);
         }
      }
   }
   else
   {
      for(int i=0;i<per;i++)
      { 
      weight += (per - i);
      sum    += array[i]*(per - i);
      }   
   }
   
   if(weight > 0) lwma = sum/weight; else lwma = 0; 
      
   return(lwma);
}

 
// MA_Method=4: SineWMA - Sine Weighted Moving Average
double SineWMAOnArray(const double& array[],int per,int bar)
{
   double swma = EMPTY_VALUE, sum = 0, weight = 0;
   
   if(bar > per)
   {
      for(int i=0;i<per;i++)
      { 
      weight += MathSin(pi*(i + 1)/(per + 1));
      sum    += array[bar-i]*MathSin(pi*(i + 1)/(per + 1)); 
      }
   
   if(weight > 0) swma = sum/weight; else swma = 0; 
   }
   
   return(swma);
}

// MA_Method=5: TriMA - Triangular Moving Average
double TriMAOnArray(const double& array[],int per,int bar)
{
   int len = (int)MathCeil((per + 1)*0.5);
   double trima = EMPTY_VALUE, sma, sum = 0;
   
   if(bar > per)
   {  
      for(int i=0;i<len;i++) 
      {
      sma  = SMAOnArray(0,array,len,bar-i);
      sum += sma;
      } 
   trima = sum/len;
   }
   
   return(trima);
}

// MA_Method=6: LSMA - Least Square Moving Average (or EPMA, Linear Regression Line)
double LSMAOnArray(const double& array[],int per,int bar)
{   
   double lsma = EMPTY_VALUE, sum = 0;
   
   if(bar > per)
   {
   for(int i=per;i>=1;i--) sum += (i - (per + 1)/3.0)*array[bar-(per-i)];
   lsma = sum*6/(per*(per + 1));
   }
   
   return(lsma);
}

// MA_Method=7: SMMA - Smoothed Moving Average
double SMMAOnArray(const double& array[],double prev,int per,int cbars,int bar)
{
   double smma = EMPTY_VALUE;
   
   if(bar == cbars) smma = SMAOnArray(0,array,per,bar); 
   else 
   if(bar >  cbars)
   {
   double sum = 0;
   for(int i=0;i<per;i++) sum += array[bar-(i+1)];
   smma = (sum - prev + array[bar])/per;
   }
   
   return(smma);
}
                                
// MA_Method=8: HMA - Hull Moving Average by Alan Hull
double HMAOnArray(const double& array[],int per,int cbars,int bar)
{
   double hma = EMPTY_VALUE, tmparray[];
   int i, len = (int)MathSqrt(per);
   
   ArrayResize(tmparray,len);
     
   if(bar == cbars) hma = array[bar]; 
   else
   if(bar >  cbars)
   {
   for(i=0;i<len;i++) tmparray[i] = 2*LWMAOnArray(0,array,per/2,bar-i) - LWMAOnArray(0,array,per,bar-i); 
   hma = LWMAOnArray(1,tmparray,len,0); 
   }  
   
   return(hma);
}

// MA_Method=9: ZeroLagEMA - Zero-Lag Exponential Moving Average
double ZeroLagEMAOnArray(const double& price[],double prev,int per,int cbars,int bar)
{
   double zema = EMPTY_VALUE, alpha = 2.0/(1 + per); 
   double lag  = 0.5*(per - 1); 
   
   if(bar == cbars) zema = price[bar];
   else 
   if(bar >  cbars) zema = alpha*(2*price[bar] - price[bar-(int)lag]) + (1 - alpha)*prev;
   
   return(zema);
}

// MA_Method=10: DEMA - Double Exponential Moving Average by Patrick Mulloy
double DEMAOnArray(int index,int num,const double price,double per,double v,int cbars,int bar)
{
   double dema = EMPTY_VALUE, alpha = 2.0/(1+per);
   
   if(bar == cbars) {dema = price; tmp[num][index][0] = dema; tmp[num+1][index][0] = dema;}
   else 
   if(bar >  cbars) 
   {
   tmp[num  ][index][0] = tmp[num  ][index][1] + alpha*(price              - tmp[num  ][index][1]); 
   tmp[num+1][index][0] = tmp[num+1][index][1] + alpha*(tmp[num][index][0] - tmp[num+1][index][1]); 
   dema                 = tmp[num  ][index][0]*(1 + v) - tmp[num+1][index][0]*v;
   }
   
   return(dema);
}

// MA_Method=11: T3 by T.Tillson
double T3_basicOnArray(int index,int num,const double price,int per,double v,int cbars,int bar)
{
   double T3 = EMPTY_VALUE, dema1, dema2;
   
   if(bar == cbars) 
   {
   T3 = price; 
   for(int k=0;k<6;k++) tmp[num+k][index][0] = T3;
   }
   else 
   if(bar  > cbars) 
   {
   dema1 = DEMAOnArray(index,num  ,price,per,v,cbars,bar); 
   dema2 = DEMAOnArray(index,num+2,dema1,per,v,cbars,bar); 
   T3    = DEMAOnArray(index,num+4,dema2,per,v,cbars,bar);
   }
   
   return(T3);
}

// MA_Method=12: ITrend - Instantaneous Trendline by J.Ehlers
double ITrendOnArray(const double& price[],double& array[],int per,int cbars,int bar)
{
   double it = EMPTY_VALUE, alpha = 2.0/(per + 1);
   
   if(bar > cbars && array[1] != EMPTY_VALUE && array[2] != EMPTY_VALUE) it = (alpha - 0.25*alpha*alpha)*price[bar] + 0.5*alpha*alpha*price[bar-1] - (alpha - 0.75*alpha*alpha)*price[bar-2] + 2*(1 - alpha)*array[1] - (1 - alpha)*(1 - alpha)*array[2];
   else 
   if(bar <= cbars && bar > 2) it = (price[bar] + 2*price[bar-1] + price[bar-2])/4; 
   
   return(it);
}

// MA_Method=13: Median - Moving Median
double MedianOnArray(const double& price[],int per,int cbars,int bar)
{
   double median = EMPTY_VALUE, array[];
   ArrayResize(array,per);
   
   if(bar >= cbars)
   {
   for(int i=0;i<per;i++) array[i] = price[bar-i];
   ArraySort(array);
   
   double num = MathRound(0.5*(per - 1)); 
   if(MathMod(per,2) > 0) median = array[(int)num]; else median = 0.5*(array[(int)num] + array[(int)num+1]);
   }
   
   return(median); 
}

// MA_Method=14: GeoMean - Geometric Mean
double GeoMeanOnArray(const double& price[],int per,int cbars,int bar)
{
   double gmean = EMPTY_VALUE;
   
   if(bar > cbars) 
   { 
   gmean = MathPow(price[bar],1.0/per); 
   for(int i=1;i<per;i++) gmean *= MathPow(price[bar-i],1.0/per); 
   }
   else 
   if(bar == cbars) gmean = SMAOnArray(0,price,per,bar);
   
   return(gmean);
}

// MA_Method=15: REMA - Regularized EMA by Chris Satchwell 
double REMAOnArray(const double price,double& array[],int per,double lambda,int cbars,int bar)
{
   double rema = EMPTY_VALUE, alpha =  2.0/(per + 1);
      
   if(bar > cbars) rema = (array[1]*(1 + 2*lambda) + alpha*(price - array[1]) - lambda*array[2])/(1 + lambda); else rema = price;
   
   return(rema);
}
// MA_Method=16: ILRS - Integral of Linear Regression Slope 
double ILRSOnArray(const double& price[],int per,int cbars,int bar)
{
   double ilrs = EMPTY_VALUE, slope, sum = per*(per - 1)*0.5, sum2 = (per - 1)*per*(2*per - 1)/6.0;
   double sum1 = 0, sumy = 0;
   
   if(bar >= cbars)
   { 
      for(int i=0;i<per;i++)
      { 
      sum1 += i*price[bar-i];
      sumy += price[bar-i];
      }
   double num1 = per*sum1 - sum*sumy;
   double num2 = sum*sum  - per*sum2;
   
   if(num2 != 0) slope = num1/num2; else slope = 0; 
   ilrs = slope + SMAOnArray(0,price,per,bar);
   }
   
   return(ilrs);
}
// MA_Method=17: IE/2 - Combination of LSMA and ILRS 
double IE2OnArray(const double& price[],int per,int cbars,int bar)
{
   double ie = EMPTY_VALUE;
   
   if(bar >= cbars) ie = 0.5*(ILRSOnArray(price,per,cbars,bar) + LSMAOnArray(price,per,bar));
      
   return(ie); 
}

// MA_Method=18: TriMAgen - Triangular Moving Average Generalized by J.Ehlers
double TriMA_genOnArray(int mode,const double& array[],int per,int cbars,int bar)
{
   int len1 = (int)MathFloor((per + 1)*0.5);
   int len2 = (int)MathCeil ((per + 1)*0.5);
   double trimagen = EMPTY_VALUE, sum = 0;
   
   if(mode == 0)
   {
      if(bar > per)
      {
      for(int i=0;i<len2;i++) sum += SMAOnArray(0,array,len1,bar-i);
      }
   }
   else
   {
   for(int i=0;i<len2;i++) sum += SMAOnArray(1,array,len1,i);
   }
   
   trimagen = sum/len2;
   
   return(trimagen);
}

// MA_Method=19: VWMA - Volume Weighted Moving Average 
double VWMAOnArray(const double& array[],int per,int rates,int bar)
{
   double vwma = EMPTY_VALUE, Sum = 0, Weight = 0;
   long volume[];
   ArraySetAsSeries(volume,true);
   
   int volumes = CopyTickVolume(_Symbol,0,rates-1-bar,per,volume);
   if(volumes < 0) return(EMPTY_VALUE);
      
   if(bar > per)
   {
      for(int i=0;i<per;i++)
      { 
      Weight += (double)volume[i];
      Sum    += (double)array[bar-i]*volume[i];
      }
   
   if(Weight > 0) vwma = Sum/Weight; else vwma = 0; 
   }
   
   return(vwma);
} 

// MA_Method=20: JSmooth - Smoothing by Mark Jurik
double JSmoothOnArray(int index,int num,const double price,int per,double pow,int cbars,int bar)
{
   double beta = 0.45*(per - 1)/(0.45*(per - 1) + 2);
	double alpha = MathPow(beta,pow);
	
	if(bar == cbars) {tmp[num+4][index][0] = price; tmp[num+0][index][0] = price; tmp[num+1][index][0] = 0; tmp[num+2][index][0] = price; tmp[num+3][index][0] = 0;}
	else 
   if(bar >  cbars) 
   {
	tmp[num+0][index][0] = (1 - alpha)*price + alpha*tmp[num+0][index][1];
	tmp[num+1][index][0] = (price - tmp[num+0][index][0])*(1 - beta) + beta*tmp[num+1][index][1];
	tmp[num+2][index][0] = tmp[num+0][index][0] + tmp[num+1][index][0];
	tmp[num+3][index][0] = (tmp[num+2][index][0] - tmp[num+4][index][1])*MathPow((1 - alpha),2) + MathPow(alpha,2)*tmp[num+3][index][1];
	tmp[num+4][index][0] = tmp[num+4][index][1] + tmp[num+3][index][0]; 
   }
  
   return(tmp[num+4][index][0]);
}

// MA_Method=21: SMA_eq     - Simplified SMA
double SMA_eqOnArray(const double& price[],double& array[],int per,int cbars,int bar)
{
   double sma = EMPTY_VALUE;
   
   if(bar == cbars) sma = SMAOnArray(0,price,per,bar);
   else 
   if(bar >  cbars) sma = (price[bar] - price[bar-per])/per + array[1]; 
     
   return(sma);
}                 
// MA_Method=22: ALMA by Arnaud Legoux / Dimitris Kouzis-Loukas / Anthony Cascino
double ALMAOnArray(const double& array[],int per,double offset,double sigma,int bar)
{
   double alma = EMPTY_VALUE, m = MathFloor(offset*(per - 1)), s = per/sigma;
	double w, sum = 0, wsum = 0;		
	
	if(bar > per)
	{
	   for(int i=0;i<per;i++) 
	   {
	   w     = MathExp(-((i - m)*(i - m))/(2*s*s));
      wsum += w;
      sum  += array[bar-(per-1-i)]*w;  
      }
   
   if(wsum != 0) alma = sum/wsum; 
   }
    
   return(alma);
}
  
// MA_Method=23: TEMA - Triple Exponential Moving Average by Patrick Mulloy
double TEMAOnArray(int index,const double price,int per,double v,int cbars,int bar)
{
   double alpha = 2.0/(per + 1);
	
	if(bar == cbars) for(int k=0;k<4;k++) tmp[k][index][0] = price; 
	else 
   if(bar >  cbars) 
   {
	tmp[0][index][0] = tmp[0][index][1] + alpha*(price            - tmp[0][index][1]);
	tmp[1][index][0] = tmp[1][index][1] + alpha*(tmp[0][index][0] - tmp[1][index][1]);
	tmp[2][index][0] = tmp[2][index][1] + alpha*(tmp[1][index][0] - tmp[2][index][1]);
	tmp[3][index][0] = tmp[0][index][0] + v*(tmp[0][index][0] + v*(tmp[0][index][0]-tmp[1][index][0]) - tmp[1][index][0] - v*(tmp[1][index][0] - tmp[2][index][0])); 
	}
   
   return(tmp[3][index][0]);
}
// MA_Method=24: T3 by T.Tillson (correct version) 
double T3OnArray(int index,int num,const double price,int per,double v,int cbars,int bar)
{
   double T3 = EMPTY_VALUE, len = MathMax((per + 5.0)/3.0 - 1,1), dema1, dema2;
   
   if(bar == cbars) 
   {
   T3 = price; 
   for(int k=0;k<6;k++) tmp[num+k][index][0] = T3;
   }
   else 
   if(bar > cbars) 
   {
   dema1 = DEMAOnArray(index,num  ,price,len,v,cbars,bar); 
   dema2 = DEMAOnArray(index,num+2,dema1,len,v,cbars,bar); 
   T3    = DEMAOnArray(index,num+4,dema2,len,v,cbars,bar);
   }
      
   return(T3);
}

// MA_Method=25: Laguerre filter by J.Ehlers
double LaguerreOnArray(int index,const double price,int per,int order,int cbars,int bar)
{
   double laguerre = EMPTY_VALUE, gamma = 1 - 10.0/(per + 9);
   double aprice[];
   
   ArrayResize(aprice,order);
   ArrayInitialize(aprice,price);
   
   if(bar >= cbars)
   {
      for(int i=0;i<order;i++)
      {
         if(bar == cbars) tmp[i][index][0] = price;
         else
         if(bar >  cbars) 
         {
            if(i == 0) tmp[i][index][0] = (1 - gamma)*price + gamma*tmp[i][index][1];
            else
            tmp[i][index][0] = -gamma*tmp[i-1][index][0] + tmp[i-1][index][1] + gamma*tmp[i][index][1];
      
         aprice[i] = tmp[i][index][0];
         }
      }
   if(bar > cbars) laguerre = TriMA_genOnArray(1,aprice,order,cbars,0); else laguerre = price; 
   }
   
   return(laguerre);
}

// MA_Method=26:  MD - McGinley Dynamic
double McGinleyOnArray(const double price,double prev,int per,int cbars,int bar)
{
   double md = EMPTY_VALUE;
   
   if(bar == cbars) md = price;
   else 
   if(bar >  cbars) 
      if(prev != 0) md = prev + (price - prev)/(per*MathPow(price/prev,4)/2); 
      else md = price;

   return(md);
}

// MA_Method=27: BF2P - Two-pole modified Butterworth filter
double BF2POnArray(const double& price[],double& array[],int per,int cbars,int bar)
{
   double bf2p = EMPTY_VALUE;
   double a  = MathExp(-1.414*pi/per);
   double b  = 2*a*MathCos(1.414*1.25*pi/per);
   double c2 = b;
   double c3 = -a*a;
   double c1 = 1 - c2 - c3;
   
   if(bar >= cbars && array[1] != EMPTY_VALUE && array[2] != EMPTY_VALUE) bf2p = c1*(price[bar] + 2*price[bar-1] + price[bar-2])/4 + c2*array[1] + c3*array[2]; 
   else if(bar < cbars && bar > 2) bf2p = (price[bar] + 2*price[bar-1] + price[bar-2])/4;
  
   return(bf2p);
}

// MA_Method=28: BF3P - Three-pole modified Butterworth filter
double BF3POnArray(const double& price[],double& array[],int per,int cbars,int bar)
{
   double bf3p = 0;
   double a  = MathExp(-pi/per);
   double b  = 2*a*MathCos(1.738*pi/per);
   double c  = a*a;
   double d2 = b + c;
   double d3 = -(c + b*c);
   double d4 = c*c;
   double d1 = 1 - d2 - d3 - d4;
   
   if(bar >= cbars && array[1] != EMPTY_VALUE && array[2] != EMPTY_VALUE && array[3] != EMPTY_VALUE) bf3p = d1*(price[bar] + 3*price[bar-1] + 3*price[bar-2] + price[bar-3])/8 + d2*array[1] + d3*array[2] + d4*array[3];
   else if(bar < cbars && bar > 2) bf3p = (price[bar] + 3*price[bar-1] + 3*price[bar-2] + price[bar-3])/8;
   
   return(bf3p);
}

// MA_Method=29: SuperSmu - SuperSmoother filter
double SuperSmuOnArray(const double& price[],double& array[],int per,int cbars,int bar)
{
   double supsm = 0;
   double a  = MathExp(-1.414*pi/per);
   double b  = 2*a*MathCos(1.414*pi/per);
   double c2 = b;
   double c3 = -a*a;
   double c1 = 1 - c2 - c3;
   
   if(bar >= cbars && array[1] != EMPTY_VALUE && array[2] != EMPTY_VALUE) supsm = c1*(price[bar] + price[bar-1])/2 + c2*array[1] + c3*array[2]; 
   else if(bar < cbars && bar > 2) supsm = (price[bar] + price[bar-1])/2;
   
   return(supsm);
}

// MA_Method=30: Decycler - Simple Decycler by J.Ehlers
double DecyclerOnArray(const double& price[],double& hp[],int per,int cbars,int bar)
{
   double alpha1 = (MathCos(1.414*pi/per) + MathSin(1.414*pi/per) - 1)/MathCos(1.414*pi/per);
  
   
   if(bar < cbars - 4) return(0); 
   
   hp[0] = (1 - alpha1/2)*(1 - alpha1/2)*(price[bar] - 2*price[bar-1] + price[bar-2]) + 2*(1 - alpha1)*hp[1] - (1 - alpha1)*(1 - alpha1)*hp[2];		
     
   return(hp[0]);
}

// MA_Method=31: eVWMA - Elastic Volume Weighted Moving Average by C.Fries
double eVWMAOnArray(const double price,double prev,int per,int cbars,int rates,int bar)
{
   double evwma = EMPTY_VALUE;
   long volume[];
   ArraySetAsSeries(volume,true);
   
   int volumes = CopyTickVolume(_Symbol,0,rates-1-bar,per,volume);
   if(volumes < 0) return(EMPTY_VALUE);
   
   if(bar == cbars) evwma = price;
   else 
   if(bar >  cbars) 
   {
   double max = 0;
   for(int i=0;i<per;i++) max = MathMax(max,volume[i]);
      
   double diff = 3*max - volume[0];
          
      if(diff < 0) evwma = prev;
      else 
      evwma = (diff*prev + volume[0]*price)/(3*max); 
   }
   
   return(evwma);
}


// MA_Method=EWMA - Exponential Weighted Moving Average 
double EWMAOnArray(const double& array[],int per,int bar)
{
   double sum = 0, weight = 0, alpha = 2.0/(1 + per);
   
   if(bar >= ma_period)
      for(int i=0;i<ma_period;i++)
      { 
      weight += alpha*MathPow(1-alpha,i);
      sum    += array[bar-i]*alpha*MathPow(1-alpha,i);
      }
      
   if(weight > 0) return(sum/weight); else return(0); 
} 

// MA_Method=DsEMA - Double Smoothed EMA
double DsEMAOnArray(int index,const double price,int per,int cbars,int bar)
{
   double alpha = 4.0/(per + 3);
   
   if(bar == cbars) 
   {
   for(int k=0;k<2;k++) tmp[k][index][0] = price;
   return(price);
   }
   else 
   if(bar > cbars) 
   {
   tmp[0][index][0] = tmp[0][index][1] + alpha*(price            - tmp[0][index][1]);
	tmp[1][index][0] = tmp[1][index][1] + alpha*(tmp[0][index][0] - tmp[1][index][1]);
   }
      
   return(tmp[1][index][0]);
}

// MA_Method=TsEMA - Triple Smoothed EMA
double TsEMAOnArray(int index,const double price,int per,int cbars,int bar)
{
   double alpha = 6.0/(per + 5);
   
   if(bar == cbars) 
   {
   for(int k=0;k<3;k++) tmp[k][index][0] = price;
   return(price);
   }
   else 
   if(bar > cbars) 
   {
   tmp[0][index][0] = tmp[0][index][1] + alpha*(price            - tmp[0][index][1]);
	tmp[1][index][0] = tmp[1][index][1] + alpha*(tmp[0][index][0] - tmp[1][index][1]);
   tmp[2][index][0] = tmp[2][index][1] + alpha*(tmp[1][index][0] - tmp[2][index][1]);
   }
      
   return(tmp[2][index][0]);
}

// MA_Method=VEMA - Volume-weighted Exponential Moving Average(V-EMA)
double VEMAOnArray(int index,double price,int per,int cbars,int rates,int bar)
{
   double vema = price, alpha = 2.0/(per + 1);
   long volume[];
   ArraySetAsSeries(volume,true);
   
   int volumes = CopyTickVolume(_Symbol,0,rates-1-bar,1,volume);
   if(volumes < 0) return(EMPTY_VALUE);
   
   if(bar == cbars) 
   {
   tmp[0][index][0] = (1 - alpha)*price*volume[0];
   tmp[1][index][0] = (1 - alpha)*volume[0];
   }
   else 
   if(bar > cbars) 
   {
   tmp[0][index][0] = tmp[0][index][1] + alpha*(price*volume[0] - tmp[0][index][1]);
	tmp[1][index][0] = tmp[1][index][1] + alpha*(volume[0]       - tmp[1][index][1]);
   }
  
   if(tmp[1][index][0] > 0) vema = tmp[0][index][0]/tmp[1][index][0];
   
   return(vema);
}


datetime prevhabar[1];
double   haClose[1][2], haOpen[1][2], haHigh[1][2], haLow[1][2];

double HeikenAshi(int index,int price,double close,double open,double high, double low,int cbars,int bar)
{ 
   if(prevhabar[index] != bar)
   {
   haClose[index][1] = haClose[index][0];
   haOpen [index][1] = haOpen [index][0];
   haHigh [index][1] = haHigh [index][0];
   haLow  [index][1] = haLow  [index][0];
   prevhabar[index]  = bar;
   }
   
   if(bar == cbars) 
   {
   haClose[index][0] = close;
   haOpen [index][0] = open;
   haHigh [index][0] = high;
   haLow  [index][0] = low ;
   }
   else
   {
   haClose[index][0] = (open + high + low + close)/4;
   haOpen [index][0] = (haOpen[index][1] + haClose[index][1])/2;
   haHigh [index][0] = MathMax(high,MathMax(haOpen[index][0],haClose[index][0]));
   haLow  [index][0] = MathMin(low ,MathMin(haOpen[index][0],haClose[index][0]));
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


double getPrice(int price,double close,double open,double high,double low,int bar)
{
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

string timeframeToString(ENUM_TIMEFRAMES timeframe)
{
   switch(timeframe)
   {
   case PERIOD_CURRENT  : return("Current");
   case PERIOD_M1       : return("M1");   
   case PERIOD_M2       : return("M2");
   case PERIOD_M3       : return("M3");
   case PERIOD_M4       : return("M4");
   case PERIOD_M5       : return("M5");      
   case PERIOD_M6       : return("M6");
   case PERIOD_M10      : return("M10");
   case PERIOD_M12      : return("M12");
   case PERIOD_M15      : return("M15");
   case PERIOD_M20      : return("M20");
   case PERIOD_M30      : return("M30");
   case PERIOD_H1       : return("H1");
   case PERIOD_H2       : return("H2");
   case PERIOD_H3       : return("H3");
   case PERIOD_H4       : return("H4");
   case PERIOD_H6       : return("H6");
   case PERIOD_H8       : return("H8");
   case PERIOD_H12      : return("H12");
   case PERIOD_D1       : return("D1");
   case PERIOD_W1       : return("W1");
   case PERIOD_MN1      : return("MN1");      
   default              : return("Current");
   }
}


datetime prevnbtime;

bool isNewBar(ENUM_TIMEFRAMES timeframe)
{
   bool res = false;
   
   if(tf >= 0)
   {
      if(iTime(NULL,timeframe,0) != prevnbtime)
      {
      res   = true;
      prevnbtime = iTime(NULL,timeframe,0);
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




