/****************************************************************************

    DRC: Digital Room Correction
    Copyright (C) 2002-2017 Denis Sbragion

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

		You can contact the author on Internet at the following address:

				d.sbragion@neomerica.it

		This program uses the FFT  routines from  Takuya Ooura and the GNU
		Scientific  Library (GSL).  Many thanks  to Takuya Ooura and the GSL
		developers for these efficient routines.

****************************************************************************/

/* Main file */

/* Inclusioni */
#include "drc.h"
#include "dsplib.h"
#include "dspwind.h"
#include "bwprefilt.h"
#include "slprefilt.h"
#include "baselib.h"
#include "level.h"
#include "cfgparse.h"
#include "convol.h"
#include "hd.h"
#include "toeplitz.h"
#include "fir.h"
#include "kirkebyfd.h"
#include "drccfg.h"
#include "cmdline.h"
#include "spline.h"
#include "psychoacoustic.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <thread>
#include <functional>
#include <fenv.h>

/* Versione corrente */
#define DRCVersion "3.2.3"
#define DRCCopyright "2002-2019"

/* Memory leaks debugger */
#ifdef DebugMLeaks
	#include "debug_new.h"
#endif

/* Header iniziale programma */
void ShowDRCHeader(void)
	{
		sputs("\nDRC " DRCVersion ": Digital Room Correction");
		sputs("Copyright (C) " DRCCopyright " Denis Sbragion");
		#ifdef UseDouble
			sputs("\nCompiled with double precision arithmetic.");
		#else
			sputs("\nCompiled with single precision arithmetic.");
		#endif
		#if defined(UseGSLFft) && defined(UseOouraFft)
			sputs("Using the GNU Scientific Library and Ooura FFT routines.");
		#else
			#ifdef UseGSLFft
				sputs("Using the GNU Scientific Library FFT routines.");
			#else
				#ifdef UseOouraFft
					sputs("Using the Ooura FFT routines.");
				#else
					sputs("Using the builtin FFT routines.");
				#endif
			#endif
		#endif
		sputs("\nThis program may be freely redistributed under the terms of");
		sputs("the GNU GPL and is provided to you as is, without any warranty");
		sputs("of any kind. Please read the file \"COPYING\" for details.");
	}

/* Descrizione uso programma */
void ShowDRCUsage(void)
	{
		sputs("Usage: DRC [--help] [Options] DRCFile");
		sputs("\nParameters:\n");
		sputs("  --help : show the full options list (long)");
		sputs("  Options: parameters overwriting options");
		sputs("  DRCFile: name of DRC configuration file\n");
		sputs("  Refer to the manual and samples for options");
		sputs("  details and file format\n");
	}

/* Main procedure */

extern void start_rest_server(int port);


#include <thread>
#include <functional>

struct DrcContext {
    DLReal * OInSig = NULL;
    DLReal * InSig;
    int PSStart;
    int PSEnd;
    DLReal * MPSig;
    DLReal * EPSig;
    int WStart1;
    int WLen1;
    int WStart2;
    int WLen2;
    int WStart3 = 0;
    int WLen3 = 0;
    DLReal * MCFilterFreqs;
    DLReal * MCFilterM;
    DLReal * MCFilterP;
    DLReal * MCFilter = NULL;
    DLReal * MCOutSig;
    int MCOutSigStart = 0;
    int MCOutSigLen;
    int MCMPFLen;
    DLReal * MPPFSig;
    int MPPFSigLen;
    DLReal * EPPFSig;
    int EPPFSigLen;
    DLReal * MPEPSig = NULL;
    int MPEPSigLen;
    DLReal * PTTConv;
    int PTTConvLen;
    int PTTConvStart;
    int PTTRefLen;
    DLReal * PTFilter;
    MKSETFType TFType = MKSETFLinearPhase;
    DLReal * ISRevSig;
    DLReal * ISMPEPSig;
    DLReal * ISRevOut = NULL;
    int ISSigLen;
    DLReal * RTSig;
    int RTSigLen;
    DLReal * PSFilterFreqs;
    DLReal * PSFilterM;
    DLReal * PSFilterP;
    DLReal * PSFilter = NULL;
    DLReal * PSOutSig;
    int PSOutSigLen;
    int PSMPFLen;
    DLReal SRMSValue;
    DLReal * TCSig;
    int TCSigLen;
    int I;
    int J;
    InterpolationType FIType = Linear;
    SLPPrefilteringType SLPType;
    BWPPrefilteringType BWPType;

    char* BCInFile = NULL;
    char* PSOutFile = NULL;
    char* PSOutFileType = NULL;
};

void process_drc(struct DrcContext& ctx) {

/* Effettua la prefinestratura del segnale */
		if (Cfg.BCPreWindowLen > 0)
			{
				sputs("Input signal prewindowing.");

				/* Verifica che la finestratura sia corretta */
				if ((Cfg.BCInitWindow / 2 - Cfg.BCPreWindowLen) < ctx.PSStart)
					sputs("!!Warning: input signal too short for correct signal prewindowing, spurious spikes may be generated.");

				for (ctx.I = 0;ctx.I < Cfg.BCInitWindow / 2 - Cfg.BCPreWindowLen;ctx.I++)
					ctx.InSig[ctx.I] = 0;
				SpacedBlackmanWindow(&ctx.InSig[Cfg.BCInitWindow / 2 - Cfg.BCPreWindowLen],Cfg.BCPreWindowLen,Cfg.BCPreWindowGap,WLeft);
			}

		/*********************************************************************************/
		/* Compensazione microfono */
		/*********************************************************************************/

		/* Verifica se abilitata */
		if (Cfg.MCFilterType[0] != 'N')
			{
				/* Verifica se si devono contare i punti filtro */
				if (Cfg.MCNumPoints == 0)
					{
						sputsp("Counting mic compensation definition file points: ",Cfg.MCPointsFile);
						Cfg.MCNumPoints = FLineCount(Cfg.MCPointsFile);
						printf("Mic compensation definition file points: %d\n",Cfg.MCNumPoints);
						fflush(stdout);
					}

				/* Alloca gli array per la generazione del filtro compensazione */
				sputs("Allocating mic compensation filter arrays.");
				ctx.MCFilterFreqs = new DLReal[Cfg.MCNumPoints];
				if (ctx.MCFilterFreqs == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}
				ctx.MCFilterM = new DLReal[Cfg.MCNumPoints];
				if (ctx.MCFilterM == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}
				ctx.MCFilterP = new DLReal[Cfg.MCNumPoints];
				if (ctx.MCFilterP == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}
				ctx.MCOutSigLen = Cfg.MCFilterLen + Cfg.BCInitWindow - 1;
				ctx.MCOutSig = new DLReal[ctx.MCOutSigLen];
				if (ctx.MCOutSig == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}

				/* Legge i punti del filtro */
				sputsp("Reading mic compensation definition file: ",Cfg.MCPointsFile);
				if (ReadPoints(Cfg.MCPointsFile,(TFMagType) Cfg.MCMagType[0],ctx.MCFilterFreqs,
					ctx.MCFilterM,ctx.MCFilterP,Cfg.MCNumPoints,Cfg.BCSampleRate) == False)
					{
						sputs("Mic compensation file input failed.");
						return;
					}

        /* Effettua l'inversione diretta */
				sputs("Mic compensation direct inversion.");
				for (ctx.I = 0;ctx.I < Cfg.MCNumPoints;ctx.I++)
          {
            ctx.MCFilterM[ctx.I] = ((DRCFloat) 1.0) / ctx.MCFilterM[ctx.I];
            ctx.MCFilterP[ctx.I] = -ctx.MCFilterP[ctx.I];
          }

				/* Verifica il tipo di interpolazione */
				switch(Cfg.MCInterpolationType[0])
					{
						case 'L':
							ctx.FIType = Linear;
						break;
						case 'G':
							ctx.FIType = Logarithmic;
						break;
						case 'R':
							ctx.FIType = SplineLinear;
						break;
						case 'S':
							ctx.FIType = SplineLogarithmic;
						break;
						case 'P':
							ctx.FIType = PCHIPLinear;
						break;
						case 'H':
							ctx.FIType = PCHIPLogarithmic;
						break;
					}

				/* Verifica il tipo di filtro da utilizzare */
				switch (Cfg.MCFilterType[0])
					{
						case 'L':
							/* Alloca gli array per il filtro */
							sputs("Allocating mic compensation filter arrays.");
							ctx.MCFilter = new DLReal[Cfg.MCFilterLen];
							if (ctx.MCFilter == NULL)
								{
									sputs("Memory allocation failed.");
									return;
								}
							for (ctx.I = 0; ctx.I < Cfg.MCFilterLen; ctx.I++)
								ctx.MCFilter[ctx.I] = 0;

							/* Calcola la dimensione richiesta per il calcolo del filtro */
							if (Cfg.MCMultExponent >= 0)
								{
									/* Calcola la potenza di due superiore a Cfg.MCFilterLen */
									for(ctx.I = 1;ctx.I <= Cfg.MCFilterLen;ctx.I <<= 1);
									ctx.I *= 1 << Cfg.MCMultExponent;
								}
							else
								ctx.I = Cfg.MCFilterLen;

							/* Calcola il filtro */
							sputs("Mic compensation FIR Filter computation...");
							if (GenericFir(ctx.MCFilter,Cfg.MCFilterLen,
								ctx.MCFilterFreqs,ctx.MCFilterM,ctx.MCFilterP,Cfg.MCNumPoints,ctx.I,ctx.FIType) == False)
								{
									sputs("FIR Filter computation failed.");
									return;
								}

							/* Effettua la finestratura del filtro */
							BlackmanWindow(ctx.MCFilter,Cfg.MCFilterLen);
						break;
						case 'M':
							/* Alloca gli array per il filtro */
							sputs("Allocating mic compensation filter arrays.");
							ctx.MCMPFLen = 1 + 2 * Cfg.MCFilterLen;
							ctx.MCFilter = new DLReal[ctx.MCMPFLen];
							if (ctx.MCFilter == NULL)
								{
									sputs("Memory allocation failed.");
									return;
								}
							for (ctx.I = 0; ctx.I < ctx.MCMPFLen; ctx.I++)
								ctx.MCFilter[ctx.I] = 0;

							/* Calcola la dimensione richiesta per il calcolo del filtro */
							if (Cfg.MCMultExponent >= 0)
								{
									/* Calcola la potenza di due superiore a Cfg.MCFilterLen */
									for(ctx.I = 1;ctx.I <= ctx.MCMPFLen;ctx.I <<= 1);
									ctx.I *= 1 << Cfg.MCMultExponent;
								}
							else
								ctx.I = ctx.MCMPFLen;

							/* Calcola il filtro */
							sputs("Mic compensation FIR Filter computation...");
							if (GenericFir(ctx.MCFilter,ctx.MCMPFLen,
								ctx.MCFilterFreqs,ctx.MCFilterM,ctx.MCFilterP,Cfg.MCNumPoints,ctx.I,ctx.FIType) == False)
								{
									sputs("FIR Filter computation failed.");
									return;
								}

							/* Alloca gli array per la deconvoluzione omomorfa */
							sputs("Allocating homomorphic deconvolution arrays.");
							ctx.MPSig = new DLReal[ctx.MCMPFLen];
							if (ctx.MPSig == NULL)
								{
									sputs("Memory allocation failed.");
									return;
								}

							/* Azzera gli array */
							for (ctx.I = 0;ctx.I < ctx.MCMPFLen;ctx.I++)
								ctx.MPSig[ctx.I] = 0;

							/* Effettua la deconvoluzione omomorfa*/
							sputs("MP mic compensation filter extraction homomorphic deconvolution stage...");
							if (CepstrumHD(ctx.MCFilter,ctx.MPSig,NULL,ctx.MCMPFLen,
								Cfg.MCMultExponent) == False)
								{
									sputs("Homomorphic deconvolution failed.");
									return;
								}

							/* Effettua la finestratura del filtro a fase minima */
							HalfBlackmanWindow(ctx.MPSig,Cfg.MCFilterLen,0,WRight);

							/* Copia il filtro a fase minima nell'array filtro */
							for (ctx.I = 0;ctx.I < Cfg.MCFilterLen;ctx.I++)
								ctx.MCFilter[ctx.I] = ctx.MPSig[ctx.I];

							/* Dealloca l'array deconvoluzione */
							delete[] ctx.MPSig;
						break;
					}

        /* Verifica se si deve salvare il filtro psicoacustico */
				if (Cfg.MCFilterFile != NULL)
					{
						/* Salva la componente MP */
						sputsp("Saving mic compensation filter: ",Cfg.MCFilterFile);
						if (WriteSignal(Cfg.MCFilterFile,ctx.MCFilter,Cfg.MCFilterLen,
							(IFileType) Cfg.MCFilterFileType[0]) == False)
							{
								sputs("Mic compensation filter save failed.");
								return;
							}
					}

				/* Convoluzione filtro segnale */
				sputs("Mic compensation FIR Filter convolution...");
				if (DFftConvolve(ctx.InSig,Cfg.BCInitWindow,ctx.MCFilter,
					Cfg.MCFilterLen,ctx.MCOutSig) == False)
					{
						perror("Convolution failed.");
						return;
					}

				/* Deallocazione array */
				delete[] ctx.MCFilter;
				delete[] ctx.InSig;

				/* Determina la dimensione della finestra di uscita */
				if (Cfg.MCOutWindow > 0)
					{
						/* Verifica il tipo di filtro */
						switch (Cfg.MCFilterType[0])
							{
								case 'L':
									/* Determina la finestratura filtro */
									ctx.MCOutSigStart = (ctx.MCOutSigLen - Cfg.MCOutWindow) / 2;
									ctx.MCOutSigLen = Cfg.MCOutWindow;

									/* Effetua la finestratura filtro */
									sputs("Mic compensated signal windowing.");
									BlackmanWindow(&ctx.MCOutSig[ctx.MCOutSigStart],ctx.MCOutSigLen);
								break;
								case 'M':
									/* Determina la finestratura filtro */
									ctx.MCOutSigStart = (Cfg.BCInitWindow - Cfg.MCOutWindow) / 2;
									ctx.MCOutSigLen = Cfg.MCOutWindow;

									/* Effetua la finestratura filtro */
									sputs("Mic compensated signal windowing.");
									BlackmanWindow(&ctx.MCOutSig[ctx.MCOutSigStart],ctx.MCOutSigLen);
								break;
							}
					}
				else
					{
						/* Verifica il tipo di filtro */
						switch (Cfg.MCFilterType[0])
							{
								case 'L':
									/* Determina la finestratura filtro */
									ctx.MCOutSigStart = 0;
									ctx.PSStart += Cfg.MCFilterLen / 2;
									ctx.PSEnd += Cfg.MCFilterLen / 2;
								break;
								case 'M':
									/* Determina la finestratura filtro */
									ctx.MCOutSigStart = 0;
									ctx.MCOutSigLen = Cfg.BCInitWindow;
								break;
							}
					}

				/* Normalizzazione segnale risultante */
				if (Cfg.MCNormFactor > 0)
					{
						sputs("Mic compensated signal normalization.");
						if (SigNormalize(&ctx.MCOutSig[ctx.MCOutSigStart],ctx.MCOutSigLen,Cfg.MCNormFactor,
							(NormType) Cfg.MCNormType[0]) == False)
							{
								sputs("Normalization failed.");
								return;
							}
					}

				/* Verifica se si deve salvare il segnale compensato */
				if (Cfg.MCOutFile != NULL)
					{
						/* Salva il segnale compensato  */
						sputsp("Saving mic compensated signal: ",Cfg.MCOutFile);
						if (WriteSignal(Cfg.MCOutFile,&ctx.MCOutSig[ctx.MCOutSigStart],ctx.MCOutSigLen,
							(IFileType) Cfg.MCOutFileType[0]) == False)
							{
								sputs("Mic compensated signal save failed.");
								return;
							}
					}

				/* Deallocazione array */
				delete[] ctx.MCFilterFreqs;
				delete[] ctx.MCFilterM;
				delete[] ctx.MCFilterP;
			}
    else
      {
        /* Imposta la lunghezza del segnale */
        ctx.MCOutSigStart = 0;
        ctx.MCOutSigLen = Cfg.BCInitWindow;
        ctx.MCOutSig = ctx.InSig;
      }

    /*********************************************************************************/
		/* Salvataggio segnale convoluzione di test */
		/*********************************************************************************/

    /* Verifica se è attiva la convoluzione di test */
    if (Cfg.TCOutFile != NULL || Cfg.PTType[0] != 'N')
      {
        /* Alloca l'array per la convoluzione di test */
        sputs("Allocating test convolution signal array.");
				ctx.OInSig = new DLReal[ctx.MCOutSigLen];
				if (ctx.OInSig == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}

        /* Copia il segnale in ingresso per la convoluzione finale */
        for (ctx.I = 0,ctx.J = ctx.MCOutSigStart;ctx.I < ctx.MCOutSigLen;ctx.I++,ctx.J++)
          ctx.OInSig[ctx.I] = ctx.MCOutSig[ctx.J];
			}

    /* Calcola il valore RMS del segnale in ingresso */
		ctx.SRMSValue = GetRMSLevel(&(ctx.MCOutSig[ctx.MCOutSigStart]),ctx.MCOutSigLen);
		if (ctx.SRMSValue > ((DLReal) 0.0))
			printf("Input signal RMS level %f (%f dB).\n",(double) ctx.SRMSValue, (double) (20 * log10((double) ctx.SRMSValue)));
		else
			printf("Input signal RMS level %f (-Inf dB).\n",(double) ctx.SRMSValue);
		fflush(stdout);

    /*********************************************************************************/
		/* Dip limiting preventivo */
		/*********************************************************************************/

		/* Verifica se si deve effettuare il dip limiting */
		if (Cfg.BCDLMinGain > 0)
			{
				switch (Cfg.BCDLType[0])
					{
						/* Fase lineare */
						case 'L':
						case 'P':
							sputs("Input signal linear phase dip limiting...");
							if (C1LPDipLimit(&ctx.MCOutSig[ctx.MCOutSigStart],ctx.MCOutSigLen,Cfg.BCDLMinGain,Cfg.BCDLStart,
								Cfg.BCSampleRate,Cfg.BCDLStartFreq,Cfg.BCDLEndFreq,Cfg.BCDLType[0] == 'P',Cfg.BCDLMultExponent) == False)
								{
									sputs("Dip limiting failed.");
									return;
								}
						break;

						/* Fase minima */
						case 'M':
						case 'W':
							sputs("Input signal minimum phase dip limiting...");
							if (C1HMPDipLimit(&ctx.MCOutSig[ctx.MCOutSigStart],ctx.MCOutSigLen,Cfg.BCDLMinGain,Cfg.BCDLStart,
								Cfg.BCSampleRate,Cfg.BCDLStartFreq,Cfg.BCDLEndFreq,Cfg.BCDLType[0] == 'W',Cfg.BCDLMultExponent) == False)
								{
									sputs("Dip limiting failed.");
									return;
								}
						break;
					}
			}

		/* Verifica se si deve effettuare rinormalizzazione */
		if (Cfg.BCNormFactor > 0)
			{
				sputs("Input signal normalization.");
				if (SigNormalize(&ctx.MCOutSig[ctx.MCOutSigStart],ctx.MCOutSigLen,Cfg.BCNormFactor,
					(NormType) Cfg.BCNormType[0]) == False)
					{
						sputs("Normalization failed.");
						return;
					}
			}

		/*********************************************************************************/
		/* Deconvoluzione omomorfa */
		/*********************************************************************************/

		/* Alloca gli array per la deconvoluzione omomorfa */
		sputs("Allocating homomorphic deconvolution arrays.");
		ctx.MPSig = new DLReal[2 * ctx.MCOutSigLen];
		if (ctx.MPSig == NULL)
			{
				sputs("Memory allocation failed.");
				return;
			}
		ctx.EPSig = new DLReal[ctx.MCOutSigLen];
		if (ctx.EPSig == NULL)
			{
				sputs("Memory allocation failed.");
				return;
			}

		/* Azzera gli array */
		for (ctx.I = 0;ctx.I < 2 * ctx.MCOutSigLen;ctx.I++)
			ctx.MPSig[ctx.I] = 0;
		for (ctx.I = 0;ctx.I < ctx.MCOutSigLen;ctx.I++)
			ctx.EPSig[ctx.I] = 0;

		/* Effettua la deconvoluzione omomorfa*/
		sputs("Homomorphic deconvolution stage...");
		if (CepstrumHD(&ctx.MCOutSig[ctx.MCOutSigStart],&ctx.MPSig[ctx.MCOutSigLen / 2 - (1 - (ctx.MCOutSigLen % 2))],ctx.EPSig,
			ctx.MCOutSigLen,Cfg.HDMultExponent) == False)
			{
				sputs("Homomorphic deconvolution failed.");
				return;
			}

		/* Verifica se si deve effettuare rinormalizzazione */
		if (Cfg.HDMPNormFactor > 0)
			{
				sputs("Minimum phase component normalization.");
				if (SigNormalize(ctx.MPSig,ctx.MCOutSigLen,Cfg.HDMPNormFactor,
					(NormType) Cfg.HDMPNormType[0]) == False)
					{
						sputs("Normalization failed.");
						return;
					}
			}
		if (Cfg.HDEPNormFactor > 0)
			{
				sputs("Excess phase component normalization.");
				if (SigNormalize(ctx.EPSig,ctx.MCOutSigLen,Cfg.HDEPNormFactor,
					(NormType) Cfg.HDEPNormType[0]) == False)
					{
						sputs("Normalization failed.");
						return;
					}
			}

		/* Verifica se si deve salvare la componente MP */
		if (Cfg.HDMPOutFile != NULL)
			{
				/* Salva la componente MP */
				sputsp("Saving minimum phase component: ",Cfg.HDMPOutFile);
				if (WriteSignal(Cfg.HDMPOutFile,ctx.MPSig,ctx.MCOutSigLen,
					(IFileType) Cfg.HDMPOutFileType[0]) == False)
					{
						sputs("Minimum phase component save failed.");
						return;
					}
			}

		/* Verifica se si deve salvare la componente EP */
		if (Cfg.HDEPOutFile != NULL)
			{
				/* Salva la componente EP */
				sputsp("Saving excess phase component: ",Cfg.HDEPOutFile);
				if (WriteSignal(Cfg.HDEPOutFile,ctx.EPSig,ctx.MCOutSigLen,
					(IFileType) Cfg.HDEPOutFileType[0]) == False)
					{
						sputs("Excess phase component save failed.");
						return;
					}
			}

		/* Dealloca il segnale di ingresso */
		delete[] ctx.MCOutSig;

		/*********************************************************************************/
		/* Prefiltratura componente MP */
		/*********************************************************************************/

		/* Alloca l'array per il segnale MP prefiltrato */
		sputs("Allocating minimum phase component prefiltering array.");
		ctx.MPPFSigLen = Cfg.MPLowerWindow + Cfg.MPFilterLen - 1;
		ctx.MPPFSig = new DLReal[ctx.MPPFSigLen];
		if (ctx.MPPFSig == NULL)
			{
				sputs("Memory allocation failed.");
				return;
			}

		/* Azzera l'array */
		for (ctx.I = 0;ctx.I < ctx.MPPFSigLen;ctx.I++)
			ctx.MPPFSig[ctx.I] = 0;

		/* Calcola il punto iniziale finestra */
		ctx.WStart1 = (ctx.MCOutSigLen - Cfg.MPLowerWindow) / 2;

		/* Verifica il tipo di funzione di prefiltratura */
		if (Cfg.MPPrefilterFctn[0] == 'P')
			{
				/* Proporzionale */
				ctx.SLPType = SLPProportional;
				ctx.BWPType = BWPProportional;
			}
		else
			{
				/* Bilineare */
				ctx.SLPType = SLPBilinear;
				ctx.BWPType = BWPBilinear;
			}

		/* Prefiltratura componente MP */
		switch (Cfg.MPPrefilterType[0])
			{
				case 'B':
					sputs("Minimum phase component band windowing.");

					/* Verifica che la finestratura sia corretta */
					if (((Cfg.BCPreWindowLen == 0) && (ctx.WStart1 < ctx.PSStart)) || ((ctx.WStart1 + Cfg.MPLowerWindow) > ctx.PSEnd))
						sputs("!!Warning: input signal too short for correct signal prefiltering, spurious spikes may be generated.");

					BWPreFilt(&ctx.MPSig[ctx.WStart1],Cfg.MPLowerWindow,Cfg.MPUpperWindow,
						Cfg.MPFilterLen,Cfg.MPBandSplit,Cfg.MPWindowExponent,
						Cfg.BCSampleRate,Cfg.MPStartFreq,Cfg.MPEndFreq,Cfg.MPWindowGap,
						ctx.MPPFSig,WFull,ctx.BWPType);
				break;

				case 'b':
					sputs("Minimum phase component single side band windowing.");

					/* Verifica che la finestratura sia corretta */
					if ((ctx.WStart1 + Cfg.MPLowerWindow) > ctx.PSEnd)
						sputs("!!Warning: input signal too short for correct signal prefiltering, spurious spikes may be generated.");

					BWPreFilt(&ctx.MPSig[ctx.WStart1],Cfg.MPLowerWindow,Cfg.MPUpperWindow,
						Cfg.MPFilterLen,Cfg.MPBandSplit,Cfg.MPWindowExponent,
						Cfg.BCSampleRate,Cfg.MPStartFreq,Cfg.MPEndFreq,Cfg.MPWindowGap,
						ctx.MPPFSig,WRight,ctx.BWPType);
				break;

				case 'S':
					sputs("Minimum phase component sliding lowpass prefiltering.");

					/* Verifica che la finestratura sia corretta */
					if (((Cfg.BCPreWindowLen == 0) && (ctx.WStart1 < ctx.PSStart)) || ((ctx.WStart1 + Cfg.MPLowerWindow) > ctx.PSEnd))
						sputs("!!Warning: input signal too short for correct signal prefiltering, spurious spikes may be generated.");

					SLPreFilt(&ctx.MPSig[ctx.WStart1],Cfg.MPLowerWindow,Cfg.MPUpperWindow,
						Cfg.MPFilterLen,Cfg.MPBandSplit,Cfg.MPWindowExponent,
						Cfg.BCSampleRate,Cfg.MPStartFreq,Cfg.MPEndFreq,Cfg.MPWindowGap,
						Cfg.MPFSharpness,ctx.MPPFSig,WFull,ctx.SLPType);
				break;

				case 's':
					sputs("Minimum phase component single side sliding lowpass prefiltering.");

					/* Verifica che la finestratura sia corretta */
					if ((ctx.WStart1 + Cfg.MPLowerWindow) > ctx.PSEnd)
						sputs("!!Warning: input signal too short for correct signal prefiltering, spurious spikes may be generated.");

					SLPreFilt(&ctx.MPSig[ctx.WStart1],Cfg.MPLowerWindow,Cfg.MPUpperWindow,
						Cfg.MPFilterLen,Cfg.MPBandSplit,Cfg.MPWindowExponent,
						Cfg.BCSampleRate,Cfg.MPStartFreq,Cfg.MPEndFreq,Cfg.MPWindowGap,
						Cfg.MPFSharpness,ctx.MPPFSig,WRight,ctx.SLPType);
				break;
			}

		/* Dealloca la componente MP */
		delete[] ctx.MPSig;

		/* Calcola la dimensione per la finestratura finale */
		if (Cfg.MPPFFinalWindow > 0)
			{
				ctx.WStart1 = (ctx.MPPFSigLen - Cfg.MPPFFinalWindow) / 2;
				ctx.WLen1 = Cfg.MPPFFinalWindow;
			}
		else
			{
				ctx.WStart1 = 0;
				ctx.WLen1 = ctx.MPPFSigLen;
			}

		/*********************************************************************************/
		/* Dip limiting */
		/*********************************************************************************/

		/* Verifica se si deve effettuare il dip limiting */
		if (Cfg.DLMinGain > 0)
			{
				switch (Cfg.DLType[0])
					{
						/* Fase lineare */
						case 'L':
						case 'P':
							sputs("MP signal linear phase dip limiting...");
							if (C1LPDipLimit(&ctx.MPPFSig[ctx.WStart1],ctx.WLen1,Cfg.DLMinGain,Cfg.DLStart,
								Cfg.BCSampleRate,Cfg.DLStartFreq,Cfg.DLEndFreq,Cfg.DLType[0] == 'P',Cfg.DLMultExponent) == False)
								{
									sputs("Dip limiting failed.");
									return;
								}
						break;

						/* Fase minima */
						case 'M':
						case 'W':
							sputs("MP signal minimum phase dip limiting...");
							if (C1HMPDipLimit(&ctx.MPPFSig[ctx.WStart1],ctx.WLen1,Cfg.DLMinGain,Cfg.DLStart,
								Cfg.BCSampleRate,Cfg.DLStartFreq,Cfg.DLEndFreq,Cfg.DLType[0] == 'W',Cfg.DLMultExponent) == False)
								{
									sputs("Dip limiting failed.");
									return;
								}
						break;
					}
			}

		/* Verifica se deve essere effettuata la rinormalizzazione MP */
		if (Cfg.MPHDRecover[0] == 'Y')
			{
				/* Alloca gli array per la deconvoluzione omomorfa */
				sputs("Allocating homomorphic deconvolution arrays.");
				ctx.MPSig = new DLReal[2 * ctx.WLen1];
				if (ctx.MPSig == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}

				/* Azzera gli array */
				for (ctx.I = 0;ctx.I < 2 * ctx.WLen1;ctx.I++)
					ctx.MPSig[ctx.I] = 0;

				/* Controlla se si deve preservare la componente EP della fase minima */
				if (Cfg.MPEPPreserve[0] == 'Y')
					{
						ctx.MPEPSig = new DLReal[ctx.WLen1];
						if (ctx.MPEPSig == NULL)
							{
								sputs("Memory allocation failed.");
								return;
							}

						/* Azzera gli array */
						for (ctx.I = 0;ctx.I < ctx.WLen1;ctx.I++)
							ctx.MPEPSig[ctx.I] = 0;
					}
				else
					ctx.MPEPSig = NULL;

				/* Effettua la deconvoluzione omomorfa*/
				sputs("MP Recover homomorphic deconvolution stage...");
				if (CepstrumHD(&ctx.MPPFSig[ctx.WStart1],&ctx.MPSig[ctx.WLen1 / 2 - (1 - (ctx.WLen1 % 2))],ctx.MPEPSig,
					ctx.WLen1,Cfg.MPHDMultExponent) == False)
					{
						sputs("Homomorphic deconvolution failed.");
						return;
					}

				/* Ricopia la componente MP nell'array originale */
				for (ctx.I = 0,ctx.J = ctx.WStart1;ctx.I < ctx.WLen1;ctx.I++,ctx.J++)
					ctx.MPPFSig[ctx.J] = ctx.MPSig[ctx.I];

				/* Dealloca l'array deconvoluzione */
				delete[] ctx.MPSig;
			}

		/* Verifica se si deve effettuare la finestratura finale */
		if (Cfg.MPPFFinalWindow > 0)
			{
				sputs("Minimum phase component final windowing.");
				BlackmanWindow(&ctx.MPPFSig[ctx.WStart1],ctx.WLen1);

				/* Controlla se si deve preservare la componente EP della fase minima */
				if (Cfg.MPHDRecover[0] == 'Y' && Cfg.MPEPPreserve[0] == 'Y')
					/* Effettua la finestratura della componente EP */
					BlackmanWindow(ctx.MPEPSig,ctx.WLen1);
			}

		/* Verifica se si deve effettuare rinormalizzazione */
		if (Cfg.MPPFNormFactor > 0)
			{
				sputs("Minimum phase component normalization.");
				if (SigNormalize(&ctx.MPPFSig[ctx.WStart1],ctx.WLen1,Cfg.MPPFNormFactor,
					(NormType) Cfg.MPPFNormType[0]) == False)
					{
						sputs("Normalization failed.");
						return;
					}
			}

		/* Verifica se si deve salvare la componente MP finestrata */
		if (Cfg.MPPFOutFile != NULL)
			{
				/* Salva la componente MP */
				sputsp("Saving minimum phase component: ",Cfg.MPPFOutFile);
				if (WriteSignal(Cfg.MPPFOutFile,&ctx.MPPFSig[ctx.WStart1],ctx.WLen1,
					(IFileType) Cfg.MPPFOutFileType[0]) == False)
					{
						sputs("Minimum phase component save failed.");
						return;
					}
			}

		/*********************************************************************************/
		/* Prefiltratura componente EP */
		/*********************************************************************************/

		/* Controlla se si deve preservare la componente EP della fase minima */
		if (Cfg.MPHDRecover[0] == 'Y' && Cfg.MPEPPreserve[0] == 'Y')
			{
				/* Alloca l'array per la convoluzione */
				sputs("Allocating minimum phase EP recovering arrays.");
				ctx.MPEPSigLen = ctx.MCOutSigLen + ctx.WLen1 - 1;
				ctx.EPPFSig = new DLReal[ctx.MPEPSigLen];
				if (ctx.EPSig == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}

				/* Effettua la convoluzione */
				sputs("Minimum phase EP recovering...");
				if (DFftConvolve(ctx.MPEPSig,ctx.WLen1,ctx.EPSig,ctx.MCOutSigLen,ctx.EPPFSig) == False)
					{
						sputs("Convolution failed.");
						return;
					}

				/* Recupera la componente EP */
				for (ctx.I = 0,ctx.J = ctx.WLen1 / 2;ctx.I < ctx.MCOutSigLen;ctx.I++, ctx.J++)
					ctx.EPSig[ctx.I] = ctx.EPPFSig[ctx.J];

				/* Dealloca l'array temporaneo convoluzione */
				delete[] ctx.MPEPSig;
				delete[] ctx.EPPFSig;
			}


		/* Alloca l'array per il segnale EP prefiltrato */
		sputs("Allocating excess phase component prefiltering array.");
		ctx.EPPFSigLen = Cfg.EPLowerWindow + Cfg.EPFilterLen - 1;
		ctx.EPPFSig = new DLReal[ctx.EPPFSigLen];
		if (ctx.EPPFSig == NULL)
			{
				sputs("Memory allocation failed.");
				return;
			}

		/* Azzera l'array */
		for (ctx.I = 0;ctx.I < ctx.EPPFSigLen;ctx.I++)
			ctx.EPPFSig[ctx.I] = 0;

		/* Calcola il punto iniziale finestra */
		ctx.WStart2 = (ctx.MCOutSigLen - Cfg.EPLowerWindow) / 2;

		/* Verifica il tipo di funzione di prefiltratura */
		if (Cfg.EPPrefilterFctn[0] == 'P')
			{
				/* Proporzionale */
				ctx.SLPType = SLPProportional;
				ctx.BWPType = BWPProportional;
			}
		else
			{
				/* Bilineare */
				ctx.SLPType = SLPBilinear;
				ctx.BWPType = BWPBilinear;
			}

		/* Prefiltratura componente EP */
		switch (Cfg.EPPrefilterType[0])
			{
				case 'B':
					sputs("Excess phase component band windowing.");

					/* Verifica che la finestratura sia corretta */
					if (((Cfg.BCPreWindowLen == 0) && (ctx.WStart2 < ctx.PSStart)) || ((ctx.WStart2 + Cfg.EPLowerWindow) > ctx.PSEnd))
						sputs("!!Warning: input signal too short for correct signal prefiltering, spurious spikes may be generated.");

					BWPreFilt(&ctx.EPSig[ctx.WStart2],Cfg.EPLowerWindow,Cfg.EPUpperWindow,
						Cfg.EPFilterLen,Cfg.EPBandSplit,Cfg.EPWindowExponent,
						Cfg.BCSampleRate,Cfg.EPStartFreq,Cfg.EPEndFreq,Cfg.EPWindowGap,
						ctx.EPPFSig,WFull,ctx.BWPType);
				break;

				case 'b':
					sputs("Excess phase component single side band windowing.");

					/* Verifica che la finestratura sia corretta */
					if ((ctx.WStart2 + Cfg.EPLowerWindow) > ctx.PSEnd)
						sputs("!!Warning: input signal too short for correct signal prefiltering, spurious spikes may be generated.");

					BWPreFilt(&ctx.EPSig[ctx.WStart2],Cfg.EPLowerWindow,Cfg.EPUpperWindow,
						Cfg.EPFilterLen,Cfg.EPBandSplit,Cfg.EPWindowExponent,
						Cfg.BCSampleRate,Cfg.EPStartFreq,Cfg.EPEndFreq,Cfg.EPWindowGap,
						ctx.EPPFSig,WRight,ctx.BWPType);
				break;

				case 'S':
					sputs("Excess phase component sliding lowpass prefiltering.");

					/* Verifica che la finestratura sia corretta */
					if (((Cfg.BCPreWindowLen == 0) && (ctx.WStart2 < ctx.PSStart)) || ((ctx.WStart2 + Cfg.EPLowerWindow) > ctx.PSEnd))
						sputs("!!Warning: input signal too short for correct signal prefiltering, spurious spikes may be generated.");

					SLPreFilt(&ctx.EPSig[ctx.WStart2],Cfg.EPLowerWindow,Cfg.EPUpperWindow,
						Cfg.EPFilterLen,Cfg.EPBandSplit,Cfg.EPWindowExponent,
						Cfg.BCSampleRate,Cfg.EPStartFreq,Cfg.EPEndFreq,Cfg.EPWindowGap,
						Cfg.EPFSharpness,ctx.EPPFSig,WFull,ctx.SLPType);
				break;

				case 's':
					sputs("Excess phase component single side sliding lowpass prefiltering.");

					/* Verifica che la finestratura sia corretta */
					if ((ctx.WStart2 + Cfg.EPLowerWindow) > ctx.PSEnd)
						sputs("!!Warning: input signal too short for correct signal prefiltering, spurious spikes may be generated.");

					SLPreFilt(&ctx.EPSig[ctx.WStart2],Cfg.EPLowerWindow,Cfg.EPUpperWindow,
						Cfg.EPFilterLen,Cfg.EPBandSplit,Cfg.EPWindowExponent,
						Cfg.BCSampleRate,Cfg.EPStartFreq,Cfg.EPEndFreq,Cfg.EPWindowGap,
						Cfg.EPFSharpness,ctx.EPPFSig,WRight,ctx.SLPType);
				break;
			}

		/* Dealloca la componente EP */
		delete[] ctx.EPSig;

		/* Determina la lunghezza della componente dopo la finestratura */
		if (Cfg.EPPFFinalWindow > 0)
			{
				ctx.WStart2 = (ctx.EPPFSigLen - Cfg.EPPFFinalWindow) / 2;
				ctx.WLen2 = Cfg.EPPFFinalWindow;
			}
		else
			{
				ctx.WStart2 = 0;
				ctx.WLen2 = ctx.EPPFSigLen;
			}

		/* Verifica se si deve effettuare riappianamento */
		if (Cfg.EPPFFlatGain > 0)
			{
				switch (Cfg.EPPFFlatType[0])
					{
						case 'L':
							sputs("Excess phase component linear phase flattening...");
							LPNormFlat(&ctx.EPPFSig[ctx.WStart2],ctx.WLen2,Cfg.EPPFFlatGain,
								Cfg.EPPFOGainFactor,Cfg.EPPFFGMultExponent);
						break;

						case 'M':
							sputs("Excess phase component minimum phase flattening...");
							CMPNormFlat(&ctx.EPPFSig[ctx.WStart2],ctx.WLen2,Cfg.EPPFFlatGain,
								Cfg.EPPFOGainFactor,Cfg.EPPFFGMultExponent);
						break;

						case 'D':
							/* Alloca gli array per la deconvoluzione omomorfa */
							sputs("Allocating homomorphic deconvolution arrays.");
							ctx.EPSig = new DLReal[ctx.WLen2];
							if (ctx.EPSig == NULL)
								{
									sputs("Memory allocation failed.");
									return;
								}

							/* Azzera gli array per la deconvoluzione omomorfa */
							for (ctx.I = 0;ctx.I < ctx.WLen2;ctx.I++)
								ctx.EPSig[ctx.I] = 0;

							/* Effettua la deconvoluzione omomorfa*/
							sputs("Excess phase component homomorphic deconvolution flattening...");
							if (CepstrumHD(&ctx.EPPFSig[ctx.WStart2],NULL,ctx.EPSig,
								ctx.WLen2,Cfg.EPPFFGMultExponent) == False)
								{
									sputs("Homomorphic deconvolution failed.");
									return;
								}

							/* Copia il risultato nell'array destinazione */
							for (ctx.I = 0,ctx.J = ctx.WStart2;ctx.I < ctx.WLen2;ctx.I++,ctx.J++)
								ctx.EPPFSig[ctx.J] = ctx.EPSig[ctx.I];

							/* Dealloca gli array per la deconvoluzione omomorfa */
							delete[] ctx.EPSig;
						break;
					}
			}

		/* Verifica se si deve effettuare la finestratura finale */
		if (Cfg.EPPFFinalWindow > 0)
			{
				sputs("Excess phase component final windowing.");
				BlackmanWindow(&ctx.EPPFSig[ctx.WStart2],ctx.WLen2);
			}

		/* Verifica se si deve effettuare rinormalizzazione */
		if (Cfg.EPPFNormFactor > 0)
			{
				sputs("Excess phase component normalization.");
				if (SigNormalize(&ctx.EPPFSig[ctx.WStart2],ctx.WLen2,Cfg.EPPFNormFactor,
					(NormType) Cfg.EPPFNormType[0]) == False)
					{
						sputs("Normalization failed.");
						return;
					}
			}

		/* Verifica se si deve salvare la componente EP finestrata */
		if (Cfg.EPPFOutFile != NULL)
			{
				/* Salva la componente MP */
				sputsp("Saving excess phase component: ",Cfg.EPPFOutFile);
				if (WriteSignal(Cfg.EPPFOutFile,&ctx.EPPFSig[ctx.WStart2],ctx.WLen2,
					(IFileType) Cfg.EPPFOutFileType[0]) == False)
					{
						sputs("Excess phase component save failed.");
						return;
					}
			}

		/*********************************************************************************/
		/* Combinazione componente MP e EP */
		/*********************************************************************************/

		/* Controlla se si deve attuare la fase PC */
		if (Cfg.ISType[0] == 'L' || Cfg.PCOutFile != NULL)
			{
				/* Alloca l'array per la convoluzione MP/EP */
				sputs("Allocating MP/EP convolution array.");
				ctx.MPEPSigLen = ctx.WLen1 + ctx.WLen2 - 1;
				ctx.MPEPSig = new DLReal[ctx.MPEPSigLen];
				if (ctx.MPEPSig == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}

				/* Convoluzione MP/EP */
				sputs("MP/EP Convolution...");
				if (DFftConvolve(&ctx.MPPFSig[ctx.WStart1],ctx.WLen1,&ctx.EPPFSig[ctx.WStart2],ctx.WLen2,ctx.MPEPSig) == False)
					{
						sputs("Convolution failed.");
						return;
					}

				/* Dealloca gli array MP/EP finestrati */
				if (Cfg.ISType[0] == 'L')
					{
						delete[] ctx.MPPFSig;
						delete[] ctx.EPPFSig;
					}

				/* Finestratura segnale risultante */
				if (Cfg.PCOutWindow > 0)
					{
						sputs("MP/EP signal windowing.");
						ctx.WStart3 = (ctx.MPEPSigLen - Cfg.PCOutWindow) / 2;
						ctx.WLen3 = Cfg.PCOutWindow;
						BlackmanWindow(&ctx.MPEPSig[ctx.WStart3],ctx.WLen3);
					}
				else
					{
						ctx.WStart3 = 0;
						ctx.WLen3 = ctx.MPEPSigLen;
					}

				/* Normalizzazione segnale risultante */
				if (Cfg.PCNormFactor > 0)
					{
						sputs("MP/EP normalization.");
						if (SigNormalize(&ctx.MPEPSig[ctx.WStart3],ctx.WLen3,Cfg.PCNormFactor,
							(NormType) Cfg.PCNormType[0]) == False)
							{
								sputs("Normalization failed.");
								return;
							}
					}

				/* Verifica se si deve salvare il segnale prefinestrato */
				if (Cfg.PCOutFile != NULL)
					{
						/* Salva la componente MP */
						sputsp("Saving MP/EP signal: ",Cfg.PCOutFile);
						if (WriteSignal(Cfg.PCOutFile,&ctx.MPEPSig[ctx.WStart3],ctx.WLen3,
							(IFileType) Cfg.PCOutFileType[0]) == False)
							{
								sputs("MP/EP signal save failed.");
								return;
							}
					}

				/* Dealloca gli array */
				if (Cfg.ISType[0] != 'L')
					delete[] ctx.MPEPSig;
			}

		/*********************************************************************************/
		/* Inversione risposta all'impulso */
		/*********************************************************************************/

		/* Verifica tipo inversione */
		switch (Cfg.ISType[0])
			{
				/* Fase lineare con matrice Toeplitz */
				case 'L':
					/* Ricalcola le finestre effettive per l'inversione */
					if (Cfg.ISOutWindow > 0)
						{
							/* Finestra di inversione predefinita */
							ctx.ISSigLen = Cfg.ISOutWindow;

							/* Verifica che il segnale in ingresso sia di lunghezza adeguata */
							if (ctx.WLen3 > ctx.ISSigLen)
								{
									/* Ricalcola la lunghezza del segnale in ingresso */
									ctx.WStart3 += (ctx.WLen3 - ctx.ISSigLen) / 2;
									ctx.WLen3 = ctx.ISSigLen;

									/* Rifiniestra il segnale per riportarlo alla lunghezza dell'inversione */
									BlackmanWindow(&ctx.MPEPSig[ctx.WStart3],ctx.ISSigLen);
								}
						}
					else
						/* Adotta la finestra precedente */
						ctx.ISSigLen = ctx.WLen3;

					/* Alloca l'array per l'inversione segnale */
					sputs("Allocating delay/reverse array.");
					ctx.ISRevSig = new DLReal[ctx.ISSigLen];
					if (ctx.ISRevSig == NULL)
						{
							sputs("Memory allocation failed.");
							return;
						}
					for (ctx.I = 0;ctx.I < ctx.ISSigLen;ctx.I++)
						ctx.ISRevSig[ctx.I] = (DLReal) 0.0;

					/* Alloca l'array per l'autocorrelazione */
					sputs("Allocating autocorrelation array.");
					ctx.ISMPEPSig = new DLReal[ctx.ISSigLen];
					if (ctx.ISMPEPSig == NULL)
						{
							sputs("Memory allocation failed.");
							return;
						}
					for (ctx.I = 0;ctx.I < ctx.ISSigLen;ctx.I++)
						ctx.ISMPEPSig[ctx.I] = (DLReal) 0.0;
					for (ctx.I = ctx.WStart3,ctx.J = (ctx.ISSigLen - ctx.WLen3) / 2;ctx.I < ctx.WStart3 + ctx.WLen3;ctx.I++,ctx.J++)
						ctx.ISMPEPSig[ctx.J] =	ctx.MPEPSig[ctx.I];

					/* Dealloca il segnale composto MP/EP */
					delete[] ctx.MPEPSig;

					/* Inversione e ritardo segnale */
					sputs("Signal delay/reverse.");
					for (ctx.I = 0,ctx.J = ctx.ISSigLen - 1;ctx.I < ctx.ISSigLen;ctx.I++,ctx.J--)
						ctx.ISRevSig[ctx.J] =	ctx.ISMPEPSig[ctx.I];

					/* Calcolo autocorrelazione e setup inversione */
					sputs("Autocorrelation computation...");
					if (AutoCorrelation(ctx.ISMPEPSig,ctx.ISSigLen) == False)
						{
							sputs("Autocorrelation computation failed.");
							return;
						}
					for (ctx.I = ctx.ISSigLen / 2; ctx.I < ctx.ISSigLen; ctx.I++)
						ctx.ISMPEPSig[ctx.I] = 0;

					/* Alloca l'array per l'inversione segnale */
					sputs("Allocating inversion array.");
					ctx.ISRevOut = new DLReal[ctx.ISSigLen];
					if (ctx.ISRevOut == NULL)
						{
							sputs("Memory allocation failed.");
							return;
						}

					/* Effettua l'inversione del segnale */
					sputs("Toeplitz least square inversion...");
					if (ToeplitzSolve(ctx.ISMPEPSig,ctx.ISRevSig,ctx.ISRevOut,ctx.ISSigLen) != 0)
						{
							sputs("Inversion failed.");
							return;
						}

					/* Dealloca gli array */
					delete[] ctx.ISRevSig;
					delete[] ctx.ISMPEPSig;

					/* Reimposta la lunghezza */
					ctx.WLen3 = ctx.ISSigLen;
				break;

				/* A fase minima con pre-echo truncation */
				case 'T':
					/* Verifica la dimensione filtro richiesta */
					if (Cfg.ISOutWindow > 0)
						ctx.WLen3 = Cfg.ISOutWindow;
					else
						ctx.WLen3 = ctx.WLen1 + ctx.WLen2 - 1;

					/* Alloca l'array per l'inversione segnale */
					sputs("Allocating inversion array.");
					ctx.ISRevOut = new DLReal[ctx.WLen3];
					if (ctx.ISRevOut == NULL)
						{
							sputs("Memory allocation failed.");
							return;
						}

					/* Verifica il tipo di funzione di prefiltratura */
					if (Cfg.ISPrefilterFctn[0] == 'P')
						/* Proporzionale */
						ctx.SLPType = SLPProportional;
					else
						/* Bilineare */
						ctx.SLPType = SLPBilinear;

					/* Inversione a fase minima selettiva */
					sputs("Pre-echo truncation fast deconvolution...");
					if (PETFDInvert(&ctx.MPPFSig[ctx.WStart1],ctx.WLen1,&ctx.EPPFSig[ctx.WStart2],ctx.WLen2,ctx.ISRevOut,ctx.WLen3,
						Cfg.ISPETType[0],Cfg.ISPELowerWindow,Cfg.ISPEUpperWindow,Cfg.ISPEStartFreq,
						Cfg.ISPEEndFreq,Cfg.ISPEFilterLen,Cfg.ISPEFSharpness,Cfg.ISPEBandSplit,
						Cfg.ISPEWindowExponent,ctx.SLPType,Cfg.ISPEOGainFactor,Cfg.BCSampleRate,
						Cfg.ISSMPMultExponent) == False)
						{
							sputs("Inversion failed.");
							return;
						}

					/* Dealloca gli array MP/EP finestrati */
					delete[] ctx.MPPFSig;
					delete[] ctx.EPPFSig;
				break;
			}

		/* Finestratura segnale risultante */
		if (Cfg.ISOutWindow > 0)
			{
				sputs("Inverted signal windowing.");
				ctx.WStart2 = (ctx.WLen3 - Cfg.ISOutWindow) / 2;
				ctx.WLen2 = Cfg.ISOutWindow;
				BlackmanWindow(&ctx.ISRevOut[ctx.WStart2],ctx.WLen2);
			}
		else
			{
				ctx.WStart2 = 0;
				ctx.WLen2 = ctx.WLen3;
				BlackmanWindow(&ctx.ISRevOut[ctx.WStart2],ctx.WLen2);
			}

		/* Normalizzazione segnale risultante */
		if (Cfg.ISNormFactor > 0)
			{
				sputs("Inverted signal normalization.");
				if (SigNormalize(&ctx.ISRevOut[ctx.WStart2],ctx.WLen2,Cfg.ISNormFactor,
					(NormType) Cfg.ISNormType[0]) == False)
					{
						sputs("Normalization failed.");
						return;
					}
			}

		/* Verifica se si deve salvare il segnale invertito */
		if (Cfg.ISOutFile != NULL)
			{
				/* Salva la componente MP */
				sputsp("Saving inverted signal: ",Cfg.ISOutFile);
				if (WriteSignal(Cfg.ISOutFile,&ctx.ISRevOut[ctx.WStart2],ctx.WLen2,
					(IFileType) Cfg.ISOutFileType[0]) == False)
					{
						sputs("Inverted signal save failed.");
						return;
					}
			}

		/*********************************************************************************/
		/* Calcolo target psicoacustico */
		/*********************************************************************************/

		/* Verifica se il target psicoacustico è abilitato */
		if (Cfg.PTType[0] != 'N')
			{
				/* Alloca l'array per la convoluzione filtro e risposta */
				sputs("Allocating psychoacoustic target reference convolution array.");
				ctx.PTTConvLen = ctx.WLen2 + ctx.MCOutSigLen - 1;
				ctx.PTTConv = new DLReal[ctx.PTTConvLen];
				if (ctx.PTTConv == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}

				/* Effettua la convoluzione tra filtro e risposta */
				sputs("Psychoacoustic target reference convolution...");
				if (DFftConvolve(ctx.OInSig,ctx.MCOutSigLen,ctx.ISRevOut,ctx.WLen2,ctx.PTTConv) == False)
					{
						sputs("Convolution failed.");
						return;
					}

				/* Effettua la finestratura della convoluzione di riferimento */
				ctx.PTTRefLen = (ctx.PTTConvLen - Cfg.PTReferenceWindow) / 2;
				for (ctx.I = 0;ctx.I < ctx.PTTRefLen;ctx.I++)
					ctx.PTTConv[ctx.I] = (DLReal) 0.0;
				BlackmanWindow(&ctx.PTTConv[ctx.PTTRefLen],Cfg.PTReferenceWindow);
				for (ctx.I = (ctx.PTTRefLen + Cfg.PTReferenceWindow);ctx.I < ctx.PTTConvLen;ctx.I++)
					ctx.PTTConv[ctx.I] = (DLReal) 0.0;

				/* Verifica se si deve effettuare il dip limiting sulla risposta target */
				if (Cfg.PTDLMinGain > 0)
					{
						switch (Cfg.PTDLType[0])
							{
								/* Fase lineare */
								case 'L':
								case 'P':
									sputs("Target reference signal linear phase dip limiting...");
									if (C1LPDipLimit(&ctx.PTTConv[ctx.PTTRefLen],Cfg.PTReferenceWindow,Cfg.PTDLMinGain,Cfg.PTDLStart,
										Cfg.BCSampleRate,Cfg.PTDLStartFreq,Cfg.PTDLEndFreq,Cfg.PTDLType[0] == 'P',Cfg.PTDLMultExponent) == False)
										{
											sputs("Dip limiting failed.");
											return;
										}
								break;

								/* Fase minima */
								case 'M':
								case 'W':
									sputs("Target reference minimum phase dip limiting...");
									if (C1HMPDipLimit(&ctx.PTTConv[ctx.PTTRefLen],Cfg.PTReferenceWindow,Cfg.PTDLMinGain,Cfg.PTDLStart,
										Cfg.BCSampleRate,Cfg.PTDLStartFreq,Cfg.PTDLEndFreq,Cfg.PTDLType[0] == 'W',Cfg.PTDLMultExponent) == False)
										{
											sputs("Dip limiting failed.");
											return;
										}
								break;
							}
					}

				/* Alloca l'array per il calcolo del filtro target */
				sputs("Allocating psychoacoustic target filter array.");
				ctx.PTFilter = new DLReal[Cfg.PTFilterLen];
				if (ctx.PTFilter == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}

				/* Imposta il tipo filtro target */
				switch (Cfg.PTType[0])
					{
						case 'L':
							ctx.TFType = MKSETFLinearPhase;
						break;

						case 'M':
							ctx.TFType = MKSETFMinimumPhase;
						break;
					}

				/* Calcola il filtro target psicoacustico */
				sputs("Computing psychoacoustic target filter...");
				if (P2MKSETargetFilter(&ctx.PTTConv[ctx.PTTRefLen],Cfg.PTReferenceWindow,Cfg.BCSampleRate,
					Cfg.PTBandWidth,Cfg.PTPeakDetectionStrength,ctx.PTFilter,ctx.TFType,
					Cfg.PTMultExponent,Cfg.PTFilterLen,Cfg.PTDLMinGain,Cfg.PTDLStart,
					Cfg.BCSampleRate,Cfg.PTDLStartFreq,Cfg.PTDLEndFreq) == False)
					{
						sputs("Psychoacoustic target filter computation failed.");
						return;
					}

				/* Dealloca l'array per la convoluzione target */
				delete[] ctx.PTTConv;

				/* Verifica se si deve salvare il filtro psicoacustico */
				if (Cfg.PTFilterFile != NULL)
					{
						/* Normalizzazione segnale risultante */
						if (Cfg.PTNormFactor > 0)
							{
								sputs("Psychoacoustic target filter normalization.");
								if (SigNormalize(ctx.PTFilter,Cfg.PTFilterLen,Cfg.PTNormFactor,
									(NormType) Cfg.PTNormType[0]) == False)
									{
										sputs("Normalization failed.");
										return;
									}
							}

						/* Salva la componente MP */
						sputsp("Saving psychoacoustic target filter: ",Cfg.PTFilterFile);
						if (WriteSignal(Cfg.PTFilterFile,ctx.PTFilter,Cfg.PTFilterLen,
							(IFileType) Cfg.PTFilterFileType[0]) == False)
							{
								sputs("Psychoacoustic target filter save failed.");
								return;
							}
					}

				/* Verifica il tipo di filtro target */
				switch (ctx.TFType)
					{
						case MKSETFLinearPhase:
							ctx.PTTConvStart = 0;
							ctx.PTTConvLen = ctx.WLen2 + Cfg.PTFilterLen  - 1;
						break;

						case MKSETFMinimumPhase:
							ctx.PTTConvStart = Cfg.PTFilterLen - 1;
							ctx.PTTConvLen = ctx.WLen2 + 2 * (Cfg.PTFilterLen - 1);
						break;
					}

				/* Alloca l'array per la convoluzione filtro e target */
				sputs("Allocating psychoacoustic target correction filter convolution array.");
				ctx.PTTConv = new DLReal[ctx.PTTConvLen];
				if (ctx.PTTConv == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}
				for (ctx.I = 0;ctx.I < ctx.PTTConvStart;ctx.I++)
					ctx.PTTConv[ctx.I] = (DLReal) 0.0;

				/* Effettua la convoluzione tra filtro e target */
				sputs("Psychoacoustic target correction filter convolution...");
				if (DFftConvolve(ctx.PTFilter,Cfg.PTFilterLen,&ctx.ISRevOut[ctx.WStart2],ctx.WLen2,&ctx.PTTConv[ctx.PTTConvStart]) == False)
					{
						sputs("Convolution failed.");
						return;
					}

				/* Dealloca il filtro target */
				delete[] ctx.PTFilter;

				/* Finestratura finale filtro risultante */
				if (Cfg.PTOutWindow > 0)
					{
						sputs("Psychoacoustic target correction filter windowing.");

						ctx.WStart2 = (ctx.PTTConvLen - Cfg.PTOutWindow) / 2;
						ctx.WLen2 = Cfg.PTOutWindow;
						BlackmanWindow(&ctx.PTTConv[ctx.WStart2],ctx.WLen2);
					}
				else
					{
						ctx.WStart2 = 0;
						ctx.WLen2 = ctx.PTTConvLen;
					}

				/* Normalizzazione segnale risultante */
				if (Cfg.PTNormFactor > 0)
					{
						sputs("Psychoacoustic target correction filter normalization.");
						if (SigNormalize(&ctx.PTTConv[ctx.WStart2],ctx.WLen2,Cfg.PTNormFactor,
							(NormType) Cfg.PTNormType[0]) == False)
							{
								sputs("Normalization failed.");
								return;
							}
					}

				/* Verifica se si deve salvare il filtro correzione psicoacustico */
				if (Cfg.PTOutFile != NULL)
					{
						/* Salva il filtro correzione psicoacustico */
						sputsp("Saving psychoacoustic target correction filter: ",Cfg.PTOutFile);
						if (WriteSignal(Cfg.PTOutFile,&ctx.PTTConv[ctx.WStart2],ctx.WLen2,
							(IFileType) Cfg.PTOutFileType[0]) == False)
							{
								sputs("Psychoacoustic target correction filter save failed.");
								return;
							}
					}

				/* Dealloca e riassegna il filtro inverso */
				delete[] ctx.ISRevOut;
				ctx.ISRevOut = ctx.PTTConv;
			}

		/*********************************************************************************/
		/* Peak limiting */
		/*********************************************************************************/

		/* Controlla se si deve effettuare il peak limiting */
		if (Cfg.PLMaxGain > 0)
			{
				switch (Cfg.PLType[0])
					{
						/* Fase lineare */
						case 'L':
						case 'P':
							sputs("Linear phase peak limiting...");
							if (C1LPPeakLimit(&ctx.ISRevOut[ctx.WStart2],ctx.WLen2,Cfg.PLMaxGain,Cfg.PLStart,
								Cfg.BCSampleRate,Cfg.PLStartFreq,Cfg.PLEndFreq,Cfg.PLType[0] == 'P',Cfg.PLMultExponent) == False)
								{
									sputs("Peak limiting failed.");
									return;
								}
						break;

						/* Fase minima */
						case 'M':
						case 'W':
							sputs("Minimum phase peak limiting...");
							if (C1HMPPeakLimit(&ctx.ISRevOut[ctx.WStart2],ctx.WLen2,Cfg.PLMaxGain,Cfg.PLStart,
								Cfg.BCSampleRate,Cfg.PLStartFreq,Cfg.PLEndFreq,Cfg.PLType[0] == 'W',Cfg.PLMultExponent) == False)
								{
									sputs("Peak limiting failed.");
									return;
								}
						break;
					}
			}

		/* Effettua la finestratura finale */
		if (Cfg.PLOutWindow > 0)
			{
				ctx.WStart2 = (ctx.WLen2 - Cfg.PLOutWindow) / 2;
				ctx.WLen2 = Cfg.PLOutWindow;
				sputs("Peak limited signal final windowing.");
				BlackmanWindow(&ctx.ISRevOut[ctx.WStart2],ctx.WLen2);
			}

		/* Normalizzazione segnale risultante */
		if (Cfg.PLNormFactor > 0)
			{
				sputs("Peak limited signal normalization.");
				if (SigNormalize(&ctx.ISRevOut[ctx.WStart2],ctx.WLen2,Cfg.PLNormFactor,
					(NormType) Cfg.PLNormType[0]) == False)
					{
						sputs("Normalization failed.");
						return;
					}
			}

		/* Verifica se si deve salvare il segnale limitato */
		if (Cfg.PLOutFile != NULL)
			{
				/* Salva il segnale limitato*/
				sputsp("Saving peak limited signal: ",Cfg.PLOutFile);
				if (WriteSignal(Cfg.PLOutFile,&ctx.ISRevOut[ctx.WStart2],ctx.WLen2,
					(IFileType) Cfg.PLOutFileType[0]) == False)
					{
						sputs("Peak limited signal save failed.");
						return;
					}
			}

		/*********************************************************************************/
		/* Troncatura ringing */
		/*********************************************************************************/

		/* Controlla se è abilitata */
		if (Cfg.RTType[0] != 'N')
			{
				/* Alloca l'array per la troncatura ringing */
				sputs("Allocating ringing truncation array.");
				ctx.RTSigLen = Cfg.RTLowerWindow + Cfg.RTFilterLen - 1;
				ctx.RTSig = new DLReal[ctx.RTSigLen];
				if (ctx.RTSig == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}

				/* Azzera l'array */
				for (ctx.I = 0;ctx.I < ctx.RTSigLen;ctx.I++)
					ctx.RTSig[ctx.I] = 0;

				/* Calcola il punto iniziale finestra */
				ctx.WStart3 = (ctx.WLen2 - Cfg.RTLowerWindow) / 2;

				/* Verifica il tipo di funzione di prefiltratura */
				if (Cfg.RTPrefilterFctn[0] == 'P')
					{
						/* Proporzionale */
						ctx.SLPType = SLPProportional;
						ctx.BWPType = BWPProportional;
					}
				else
					{
						/* Bilineare */
						ctx.SLPType = SLPBilinear;
						ctx.BWPType = BWPBilinear;
					}

				/* Prefiltratura componente EP */
				switch (Cfg.RTType[0])
					{
						case 'B':
							sputs("Ringing truncation band windowing.");
							BWPreFilt(&ctx.ISRevOut[ctx.WStart3],Cfg.RTLowerWindow,Cfg.RTUpperWindow,
								Cfg.RTFilterLen,Cfg.RTBandSplit,Cfg.RTWindowExponent,
								Cfg.BCSampleRate,Cfg.RTStartFreq,Cfg.RTEndFreq,Cfg.RTWindowGap,
								ctx.RTSig,WFull,ctx.BWPType);
						break;

						case 'b':
							sputs("Ringing truncation single side band windowing.");
							BWPreFilt(&ctx.ISRevOut[ctx.WStart3],Cfg.RTLowerWindow,Cfg.RTUpperWindow,
								Cfg.RTFilterLen,Cfg.RTBandSplit,Cfg.RTWindowExponent,
								Cfg.BCSampleRate,Cfg.RTStartFreq,Cfg.RTEndFreq,Cfg.RTWindowGap,
								ctx.RTSig,WRight,ctx.BWPType);
						break;

						case 'S':
							sputs("Ringing truncation sliding lowpass filtering.");
							SLPreFilt(&ctx.ISRevOut[ctx.WStart3],Cfg.RTLowerWindow,Cfg.RTUpperWindow,
								Cfg.RTFilterLen,Cfg.RTBandSplit,Cfg.RTWindowExponent,
								Cfg.BCSampleRate,Cfg.RTStartFreq,Cfg.RTEndFreq,Cfg.RTWindowGap,
								Cfg.RTFSharpness,ctx.RTSig,WFull,ctx.SLPType);
						break;

						case 's':
							sputs("Ringing truncation single side sliding lowpass filtering.");
							SLPreFilt(&ctx.ISRevOut[ctx.WStart3],Cfg.RTLowerWindow,Cfg.RTUpperWindow,
								Cfg.RTFilterLen,Cfg.RTBandSplit,Cfg.RTWindowExponent,
								Cfg.BCSampleRate,Cfg.RTStartFreq,Cfg.RTEndFreq,Cfg.RTWindowGap,
								Cfg.RTFSharpness,ctx.RTSig,WRight,ctx.SLPType);
						break;
					}

				/* Dealloca il segnale invertito */
				delete[] ctx.ISRevOut;

				/* Determina la lunghezza della componente dopo la finestratura */
				if (Cfg.RTOutWindow > 0)
					{
						ctx.WStart2 = (ctx.RTSigLen - Cfg.RTOutWindow) / 2;
						ctx.WLen2 = Cfg.RTOutWindow;
					}
				else
					{
						ctx.WStart2 = 0;
						ctx.WLen2 = ctx.RTSigLen;
					}

				/* Verifica se si deve effettuare la finestratura finale */
				if (Cfg.RTOutWindow > 0)
					{
						sputs("Ringing truncation final windowing.");
						BlackmanWindow(&ctx.RTSig[ctx.WStart2],ctx.WLen2);
					}

				/* Verifica se si deve effettuare rinormalizzazione */
				if (Cfg.RTNormFactor > 0)
					{
						sputs("Ringing truncation normalization.");
						if (SigNormalize(&ctx.RTSig[ctx.WStart2],ctx.WLen2,Cfg.RTNormFactor,
							(NormType) Cfg.RTNormType[0]) == False)
							{
								sputs("Normalization failed.");
								return;
							}
					}

				/* Verifica se si deve salvare la troncatura ringing */
				if (Cfg.RTOutFile != NULL)
					{
						/* Salva la componente MP */
						sputsp("Saving ringing truncation: ",Cfg.RTOutFile);
						if (WriteSignal(Cfg.RTOutFile,&ctx.RTSig[ctx.WStart2],ctx.WLen2,
							(IFileType) Cfg.RTOutFileType[0]) == False)
							{
								sputs("Ringing truncation save failed.");
								return;
							}
					}

				/* Reimposta il segnale invertito */
				ctx.ISRevOut = ctx.RTSig;
			}


		/*********************************************************************************/
		/* Applicazione risposta target */
		/*********************************************************************************/

		/* Verifica se si devono contare i punti filtro */
		if (Cfg.PSNumPoints == 0)
			{
				sputsp("Counting target response definition file points: ",Cfg.PSPointsFile);
				Cfg.PSNumPoints = FLineCount(Cfg.PSPointsFile);
				printf("Target response definition file points: %d\n",Cfg.PSNumPoints);
				fflush(stdout);
			}

		/* Alloca gli array per la generazione della risposta target */
		sputs("Allocating target response arrays.");
		ctx.PSFilterFreqs = new DLReal[Cfg.PSNumPoints];
		if (ctx.PSFilterFreqs == NULL)
			{
				sputs("Memory allocation failed.");
				return;
			}
		ctx.PSFilterM = new DLReal[Cfg.PSNumPoints];
		if (ctx.PSFilterM == NULL)
			{
				sputs("Memory allocation failed.");
				return;
			}
		ctx.PSFilterP = new DLReal[Cfg.PSNumPoints];
		if (ctx.PSFilterP == NULL)
			{
				sputs("Memory allocation failed.");
				return;
			}
		ctx.PSOutSigLen = Cfg.PSFilterLen + ctx.WLen2 - 1;
		ctx.PSOutSig = new DLReal[ctx.PSOutSigLen];
		if (ctx.PSOutSig == NULL)
			{
				sputs("Memory allocation failed.");
				return;
			}

		/* Legge i punti del filtro */
		sputsp("Reading target response definition file: ",Cfg.PSPointsFile);
		if (ReadPoints(Cfg.PSPointsFile,(TFMagType) Cfg.PSMagType[0],ctx.PSFilterFreqs,
			ctx.PSFilterM,ctx.PSFilterP,Cfg.PSNumPoints,Cfg.BCSampleRate) == False)
			{
				sputs("Target response point file input failed.");
				return;
			}

		/* Verifica il tipo di interpolazione */
		switch(Cfg.PSInterpolationType[0])
			{
				case 'L':
					ctx.FIType = Linear;
				break;
				case 'G':
					ctx.FIType = Logarithmic;
				break;
				case 'R':
					ctx.FIType = SplineLinear;
				break;
				case 'S':
					ctx.FIType = SplineLogarithmic;
				break;
				case 'P':
					ctx.FIType = PCHIPLinear;
				break;
				case 'H':
					ctx.FIType = PCHIPLogarithmic;
				break;
			}

		/* Verifica il tipo di filtro da utilizzare */
		switch (Cfg.PSFilterType[0])
			{
				case 'L':
					/* Alloca gli array per il filtro */
					sputs("Allocating target filter arrays.");
					ctx.PSFilter = new DLReal[Cfg.PSFilterLen];
					if (ctx.PSFilter == NULL)
						{
							sputs("Memory allocation failed.");
							return;
						}
					for (ctx.I = 0; ctx.I < Cfg.PSFilterLen; ctx.I++)
						ctx.PSFilter[ctx.I] = 0;

					/* Calcola la dimensione richiesta per il calcolo del filtro */
					if (Cfg.PSMultExponent >= 0)
						{
							/* Calcola la potenza di due superiore a Cfg.PSFilterLen */
							for(ctx.I = 1;ctx.I <= Cfg.PSFilterLen;ctx.I <<= 1);
							ctx.I *= 1 << Cfg.PSMultExponent;
						}
					else
						ctx.I = Cfg.PSFilterLen;

					/* Calcola il filtro */
					sputs("FIR Filter computation...");
					if (GenericFir(ctx.PSFilter,Cfg.PSFilterLen,
						ctx.PSFilterFreqs,ctx.PSFilterM,ctx.PSFilterP,Cfg.PSNumPoints,ctx.I,ctx.FIType) == False)
						{
							sputs("FIR Filter computation failed.");
							return;
						}

					/* Effettua la finestratura del filtro */
					BlackmanWindow(ctx.PSFilter,Cfg.PSFilterLen);
				break;
				case 'M':
				case 'T':
					/* Alloca gli array per il filtro */
					sputs("Allocating target filter arrays.");
					ctx.PSMPFLen = 1 + 2 * Cfg.PSFilterLen;
					ctx.PSFilter = new DLReal[ctx.PSMPFLen];
					if (ctx.PSFilter == NULL)
						{
							sputs("Memory allocation failed.");
							return;
						}
					for (ctx.I = 0; ctx.I < ctx.PSMPFLen; ctx.I++)
						ctx.PSFilter[ctx.I] = 0;

					/* Calcola la dimensione richiesta per il calcolo del filtro */
					if (Cfg.PSMultExponent >= 0)
						{
							/* Calcola la potenza di due superiore a Cfg.PSFilterLen */
							for(ctx.I = 1;ctx.I <= ctx.PSMPFLen;ctx.I <<= 1);
							ctx.I *= 1 << Cfg.PSMultExponent;
						}
					else
						ctx.I = ctx.PSMPFLen;

					/* Calcola il filtro */
					sputs("FIR Filter computation...");
					if (GenericFir(ctx.PSFilter,ctx.PSMPFLen,
						ctx.PSFilterFreqs,ctx.PSFilterM,ctx.PSFilterP,Cfg.PSNumPoints,ctx.I,ctx.FIType) == False)
						{
							sputs("FIR Filter computation failed.");
							return;
						}

					/* Alloca gli array per la deconvoluzione omomorfa */
					sputs("Allocating homomorphic deconvolution arrays.");
					ctx.MPSig = new DLReal[ctx.PSMPFLen];
					if (ctx.MPSig == NULL)
						{
							sputs("Memory allocation failed.");
							return;
						}

					/* Azzera gli array */
					for (ctx.I = 0;ctx.I < ctx.PSMPFLen;ctx.I++)
						ctx.MPSig[ctx.I] = 0;

					/* Effettua la deconvoluzione omomorfa*/
					sputs("MP target response extraction homomorphic deconvolution stage...");
					if (CepstrumHD(ctx.PSFilter,ctx.MPSig,NULL,ctx.PSMPFLen,
						Cfg.PSMultExponent) == False)
						{
							sputs("Homomorphic deconvolution failed.");
							return;
						}

					/* Effettua la finestratura del filtro a fase minima */
					HalfBlackmanWindow(ctx.MPSig,Cfg.PSFilterLen,0,WRight);

					/* Copia il filtro a fase minima nell'array filtro */
					for (ctx.I = 0;ctx.I < Cfg.PSFilterLen;ctx.I++)
						ctx.PSFilter[ctx.I] = ctx.MPSig[ctx.I];

					/* Dealloca l'array deconvoluzione */
					delete[] ctx.MPSig;
				break;
			}

		/* Convoluzione filtro segnale */
		sputs("Target response FIR Filter convolution...");
		if (DFftConvolve(&ctx.ISRevOut[ctx.WStart2],ctx.WLen2,ctx.PSFilter,
			Cfg.PSFilterLen,ctx.PSOutSig) == False)
			{
				perror("Convolution failed.");
				return;
			}

		/* Deallocazione array */
		delete[] ctx.ISRevOut;
		delete[] ctx.PSFilter;

		/* Determina la dimensione della finestra di uscita */
		if (Cfg.PSOutWindow > 0)
			{
				/* Alloca l'array temporaneo per il filtro */
				ctx.PSFilter = new DLReal[Cfg.PSOutWindow];
				if (ctx.PSFilter == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}

				/* Verifica il tipo di filtro */
				switch (Cfg.PSFilterType[0])
					{
						case 'L':
							/* Determina la finestratura filtro */
							ctx.WStart2 = (ctx.PSOutSigLen - Cfg.PSOutWindow) / 2;
							ctx.WLen2 = Cfg.PSOutWindow;
							ctx.WLen3 = ctx.PSOutSigLen;

							/* Salva il filtro per la convoluzione test */
							for (ctx.I = 0,ctx.J = ctx.WStart2;ctx.I < ctx.WLen2;ctx.I++,ctx.J++)
								ctx.PSFilter[ctx.I] = ctx.PSOutSig[ctx.J];

							/* Effetua la finestratura filtro */
							sputs("Target response signal windowing.");
							BlackmanWindow(ctx.PSFilter,ctx.WLen2);
						break;
						case 'M':
							/* Determina la finestratura filtro */
							ctx.WStart2 = (ctx.WLen2 - Cfg.PSOutWindow) / 2;
							ctx.WLen3 = ctx.WLen2;
							ctx.WLen2 = Cfg.PSOutWindow;

							/* Salva il filtro per la convoluzione test */
							for (ctx.I = 0,ctx.J = ctx.WStart2;ctx.I < ctx.WLen2;ctx.I++,ctx.J++)
								ctx.PSFilter[ctx.I] = ctx.PSOutSig[ctx.J];

							/* Effetua la finestratura filtro */
							sputs("Target response signal windowing.");
							BlackmanWindow(ctx.PSFilter,ctx.WLen2);
						break;
						case 'T':
							/* Determina la finestratura filtro */
							ctx.WStart2 = (ctx.WLen2 / 2) - Cfg.ISPELowerWindow;
							ctx.WLen3 = ctx.WLen2;
							ctx.WLen2 = Cfg.PSOutWindow;

							/* Salva il filtro per la convoluzione test */
							for (ctx.I = 0,ctx.J = ctx.WStart2;ctx.I < ctx.WLen2;ctx.I++,ctx.J++)
								ctx.PSFilter[ctx.I] = ctx.PSOutSig[ctx.J];

							/* Effetua la finestratura filtro */
							sputs("Target response signal windowing.");
							HalfBlackmanWindow(ctx.PSFilter,ctx.WLen2,Cfg.ISPELowerWindow,WRight);
						break;
					}
			}
		else
			{
				/* Verifica il tipo di filtro */
				switch (Cfg.PSFilterType[0])
					{
						case 'L':
							/* Determina la finestratura filtro */
							ctx.WStart2 = 0;
							ctx.WLen2 = ctx.PSOutSigLen;
							ctx.WLen3 = ctx.PSOutSigLen;
						case 'M':
							/* Determina la finestratura filtro */
							ctx.WStart2 = 0;
							ctx.WLen3 = ctx.WLen2;
						break;
						case 'T':
							/* Determina la finestratura filtro */
							ctx.WStart2 = (ctx.WLen2 / 2) - Cfg.ISPELowerWindow;
							ctx.WLen3 = ctx.WLen2;
							ctx.WLen2 = ctx.PSOutSigLen - ctx.WStart2;
						break;
					}

				/* Alloca l'array temporaneo per il filtro */
				ctx.PSFilter = new DLReal[ctx.WLen2];
				if (ctx.PSFilter == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}

				/* Salva il filtro per la convoluzione test */
				for (ctx.I = 0,ctx.J = ctx.WStart2;ctx.I < ctx.WLen2;ctx.I++,ctx.J++)
					ctx.PSFilter[ctx.I] = ctx.PSOutSig[ctx.J];
			}

		/* Normalizzazione segnale risultante */
		if (Cfg.PSNormFactor > 0)
			{
				sputs("Target response signal normalization.");
				if (SigNormalize(ctx.PSFilter,ctx.WLen2,Cfg.PSNormFactor,
					(NormType) Cfg.PSNormType[0]) == False)
					{
						sputs("Normalization failed.");
						return;
					}
			}

		/* Verifica se si deve salvare il segnale risposta target */
		if (ctx.PSOutFile != NULL)
			{
				/* Salva la componente MP */
				sputsp("Saving Target response signal: ",ctx.PSOutFile);
				if (WriteSignal(ctx.PSOutFile,ctx.PSFilter,ctx.WLen2,
					(IFileType) ctx.PSOutFileType[0]) == False)
					{
						sputs("Target response signal save failed.");
						return;
					}
			}

		/* Deallocazione array */
		delete[] ctx.PSFilterFreqs;
		delete[] ctx.PSFilterM;
		delete[] ctx.PSFilterP;

    /* Reimposta la lunghezza filtro */
    ctx.PSOutSigLen = ctx.WLen2;

		/*********************************************************************************/
		/* Estrazione filtro a fase minima */
		/*********************************************************************************/

		/* Verifica se deve essere estratto il filtro a fase minima */
		if (Cfg.MSOutFile != NULL)
			{
				/* Alloca gli array per la deconvoluzione omomorfa */
				sputs("Allocating homomorphic deconvolution arrays.");
				ctx.PSMPFLen = Cfg.MSFilterDelay + ctx.WLen2;
				ctx.MPSig = new DLReal[ctx.PSMPFLen];
				if (ctx.MPSig == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}

				/* Azzera gli array */
				for (ctx.I = 0;ctx.I < ctx.PSMPFLen;ctx.I++)
					ctx.MPSig[ctx.I] = 0;

				/* Effettua la deconvoluzione omomorfa*/
				sputs("MP filter extraction homomorphic deconvolution stage...");
				if (CepstrumHD(&ctx.PSOutSig[ctx.WStart2],&ctx.MPSig[Cfg.MSFilterDelay],NULL,
					ctx.WLen2,Cfg.MSMultExponent) == False)
					{
						sputs("Homomorphic deconvolution failed.");
						return;
					}

				/* Verifica se si deve finestrare il filtro */
				sputs("MP filter extraction windowing.");
				if (Cfg.MSOutWindow > 0)
					{
						HalfBlackmanWindow(&ctx.MPSig[Cfg.MSFilterDelay],Cfg.MSOutWindow - Cfg.MSFilterDelay,0,WRight);
						ctx.WLen1 = Cfg.MSOutWindow;
					}
				else
					ctx.WLen1 = ctx.PSMPFLen;

				/* Normalizzazione segnale risultante */
				if (Cfg.MSNormFactor > 0)
					{
						sputs("Minimum phase filter normalization.");
						if (SigNormalize(ctx.MPSig,ctx.WLen1,Cfg.MSNormFactor,
							(NormType) Cfg.MSNormType[0]) == False)
							{
								sputs("Normalization failed.");
								return;
							}
					}

				/* Salva il il filtro a fase minima */
				sputsp("Saving MP filter signal: ",Cfg.MSOutFile);
				if (WriteSignal(Cfg.MSOutFile,ctx.MPSig,ctx.WLen1,
					(IFileType) Cfg.MSOutFileType[0]) == False)
					{
						sputs("MP filter signal save failed.");
						return;
					}

				/* Dealloca l'array deconvoluzione */
				delete[] ctx.MPSig;
			}

		/* Deallocazione array */
		delete[] ctx.PSOutSig;

		/*********************************************************************************/
		/* Convoluzione di test */
		/*********************************************************************************/

		/* Verifica se va effettuata la convoluzione finale */
		if (Cfg.TCOutFile != NULL)
			{
				/* Alloca l'array per la convoluzione finale */
				sputs("Allocating test convolution arrays.");
				ctx.TCSigLen = ctx.MCOutSigLen + ctx.PSOutSigLen - 1;
				ctx.TCSig = new DLReal[ctx.TCSigLen];
				if (ctx.TCSig == NULL)
					{
						sputs("Memory allocation failed.");
						return;
					}

				/* Effettua la convoluzione */
				sputs("Convolving input signal with target response signal...");
				if (DFftConvolve(ctx.OInSig,ctx.MCOutSigLen,ctx.PSFilter,ctx.PSOutSigLen,ctx.TCSig) == False)
					{
						sputs("Convolution failed.");
						return;
					}

				/* Calcola il valore RMS del segnale dopo la filtratura */
				ctx.SRMSValue = GetRMSLevel(ctx.TCSig,ctx.TCSigLen);
				if (ctx.SRMSValue >= 0)
					printf("Filtered signal RMS level %f (%f dB).\n",(double) ctx.SRMSValue, (double) (20 * log10((double) ctx.SRMSValue)));
				else
					printf("Filtered signal RMS level %f (-inf dB).\n",(double) ctx.SRMSValue);
				fflush(stdout);

				/* Normalizzazione segnale risultante */
				if (Cfg.TCNormFactor > 0)
					{
						sputs("Test convolution signal normalization.");
						if (SigNormalize(ctx.TCSig,ctx.TCSigLen,Cfg.TCNormFactor,
							(NormType) Cfg.TCNormType[0]) == False)
							{
								sputs("Normalization failed.");
								return;
							}
					}

				/* Calcola la dimensione in uscita */
				if (Cfg.PSFilterType[0] == 'T')
					ctx.WLen3 = ctx.MCOutSigLen + 2 * Cfg.ISPELowerWindow;
				else
					ctx.WLen3 = ctx.TCSigLen;

				/* Salva il segnale convoluzione test */
				sputsp("Saving test convolution signal: ",Cfg.TCOutFile);
				if (WriteSignal(Cfg.TCOutFile,ctx.TCSig,ctx.WLen3,
					(IFileType) Cfg.TCOutFileType[0]) == False)
					{
						sputs("Test convolution save failed.");
						return;
					}

				/* Effettua la sovrascrittura del segnale convoluzione test */
				if (Cfg.TCOWFile != NULL)
					{
						sputsp("Saving test convolution overwrite: ",Cfg.TCOWFile);

						/* Normalizzazione segnale risultante */
						if (Cfg.TCOWNormFactor > 0)
							{
								sputs("Test convolution overwrite signal normalization.");
								if (SigNormalize(ctx.TCSig,ctx.TCSigLen,Cfg.TCOWNormFactor,
									(NormType) Cfg.TCOWNormType[0]) == False)
									{
										sputs("Normalization failed.");
										return;
									}
							}

						/* Controlla il tipo di filtro */
						if (Cfg.PSFilterType[0] == 'T')
							{

								if (((ctx.MCOutSigLen / 2 + Cfg.ISPELowerWindow) - Cfg.TCOWPrewindow) < (ctx.TCSigLen - Cfg.TCOWLength))
									ctx.WLen3 = Cfg.TCOWLength;
								else
									ctx.WLen3 = ctx.TCSigLen - ((ctx.MCOutSigLen / 2 + Cfg.ISPELowerWindow) - Cfg.TCOWPrewindow);

								if (OverwriteSignal(Cfg.TCOWFile,&ctx.TCSig[(ctx.MCOutSigLen / 2 + Cfg.ISPELowerWindow) - Cfg.TCOWPrewindow],
									ctx.WLen3,Cfg.TCOWSkip,(IFileType) Cfg.TCOWFileType[0]) == False)
									{
										sputs("Test convolution overwrite failed.");
										return;
									}
							}
						else
							{
								if ((ctx.TCSigLen / 2 - Cfg.TCOWPrewindow) < (ctx.TCSigLen - Cfg.TCOWLength))
									ctx.WLen3 = Cfg.TCOWLength;
								else
									ctx.WLen3 = ctx.TCSigLen / 2 + Cfg.TCOWPrewindow;

								if (OverwriteSignal(Cfg.TCOWFile,&ctx.TCSig[ctx.TCSigLen / 2 - Cfg.TCOWPrewindow],
									ctx.WLen3,Cfg.TCOWSkip,(IFileType) Cfg.TCOWFileType[0]) == False)
									{
										sputs("Test convolution overwrite failed.");
										return;
									}
							}
					}

				/* Dealloca gli array temporanei */
				delete[] ctx.TCSig;
			}

		/* Dealloca gli array temporanei */
		if (Cfg.TCOutFile != NULL || Cfg.PTType[0] != 'N')
			delete[] ctx.OInSig;

		/* Dealloca il filtro convoluzione test */
		delete[] ctx.PSFilter;

		/* Libera la memoria della struttura di configurazione */
		CfgFree(CfgParmsDef);
		//free(DRCFile);

		/* Esecuzione completata */
		sputs("Execution completed.");

		/* Segnala la durata */
		printf("Total computing time: %lu s\n",(unsigned long int) (0 /*time(NULL) - CStart*/));
		fflush(stdout);

		
}


DrcContext ctx1;
DrcContext ctx2;

void trigger_recalculation() {
    std::thread t1(process_drc, std::ref(ctx1));
    std::thread t2;
    if (Cfg.BCInFile2 != NULL && Cfg.PSOutFile2 != NULL) {
        t2 = std::thread(process_drc, std::ref(ctx2));
    }
    t1.join();
    if (t2.joinable()) {
        t2.join();
    }
}
int main(int argc, char * argv[]) {

/* Gestione parametri recuperati dalla linea di comando */
		CmdLineType * OptData;
		char * DRCFile;

		/* Salvataggio istante di avvio */
		time_t CStart = (time_t) 0;

		/* I386 Debug only, enables all floating point exceptions traps */
		/* int em = 0x372;
		__asm__ ("fldcw %0" : : "m" (em)); */

		/* Messaggio iniziale */
		ShowDRCHeader();

		/* Empty line */
		sputs("");

		/* Controllo presenza argomenti */
		if (argc < 2)
			{
				ShowDRCUsage();
				return 0;
			}

		/* Salvataggio istante di avvio */
		CStart = time(NULL);

		/* Registra le informazioni command line sulla base della struttura di
		configurazione */
		OptData = RegisterCmdLine(CfgParmsDef);
		if (OptData == NULL)
			{
				sputs("Memory allocation failed.");
				return 1;
			}

		/* Recupera i parametri della command line */
		if (GetCmdLine(argc,argv,CfgParmsDef,OptData,&DRCFile) != 0)
			{
				sputs("\nCommand line parsing error.");
				return 1;
			}

		/* Verifica se è stato richiesto l'help */
		if (OptData->ParmSet[OptData->OptCount] == True)
			{
				/* Visualizza le opzioni disponibili a linea di comando */
				ShowDRCUsage();
				sputs("Available options:\n");
				ShowCmdLine(CfgParmsDef);

				/* Dealloca le informazioni parsing command line */
				FreeCmdLine(OptData, CfgParmsDef);

				return 0;
			}

		/* Verifica che il nome del file sia presente */
		if (DRCFile == NULL)
			{
				ShowDRCUsage();
				return 1;
			}

		/* Segnala l'avvio della procedura */
		sputsp("Input configuration file: ",DRCFile);

		/* Recupera la configurazione */
		sputs("Parsing configuration file...");
		if (CfgParse(DRCFile,CfgParmsDef,CfgSimple) <= 0)
			{
				/* Dealloca le informazioni parsing command line */
				FreeCmdLine(OptData, CfgParmsDef);
				CfgFree(CfgParmsDef);

				sputs(CfgGetLastErrorDsc());
				sputs("Configuration file parsing error.");
				return 1;
			}
		sputs("Parsing completed.");

		/* Sovrascrive la configurazione base con i parametri
		a linea di comando. */
		sputs("Adding command line options...");
		CopyCmdLineParms(OptData,CfgParmsDef);

		/* Imposta la directory base recupero file */
		if (SetupDRCCfgBaseDir(&Cfg,CfgParmsDef,OptData) > 0)
			{
				/* Dealloca le informazioni parsing command line */
				FreeCmdLine(OptData, CfgParmsDef);
				CfgFree(CfgParmsDef);

				sputs("Base configuration setup error.");
				return 1;
			}

		/* Dealloca le informazioni parsing command line */
		FreeCmdLine(OptData, CfgParmsDef);

		/* Controllo validità parametri */
		sputs("Configuration parameters check.");
		if (CheckDRCCfg(&Cfg) != 0)
			{
				/* Libera la memoria della struttura di configurazione */
				CfgFree(CfgParmsDef);
				//free(DRCFile);
				return 1;
			}

		/* Controlla se è stata definita un directory base */
		if (Cfg.BCBaseDir != NULL)
			if (strlen(Cfg.BCBaseDir) > 0)
				sputsp("Base directory: ",Cfg.BCBaseDir);

		/*********************************************************************************/
		/* Importazione iniziale risposta all'impulso */
		/*********************************************************************************/

		
        // Context 1 setup
		if (Cfg.BCImpulseCenterMode[0] == 'A')
			{
				Cfg.BCImpulseCenter = FindMaxPcm(Cfg.BCInFile,(IFileType) Cfg.BCInFileType[0]);
				if (Cfg.BCImpulseCenter < 0) return 1;
			}
		ctx1.InSig = new DLReal[Cfg.BCInitWindow];
		if (ctx1.InSig == NULL) return 1;
		if (ReadSignal(Cfg.BCInFile,ctx1.InSig,Cfg.BCInitWindow,Cfg.BCImpulseCenter,(IFileType) Cfg.BCInFileType[0],&ctx1.PSStart,&ctx1.PSEnd) == False) return 1;
        ctx1.BCInFile = Cfg.BCInFile;
        ctx1.PSOutFile = Cfg.PSOutFile;
        ctx1.PSOutFileType = Cfg.PSOutFileType;

        // Context 2 setup
        if (Cfg.BCInFile2 != NULL && Cfg.PSOutFile2 != NULL) {
            if (Cfg.BCImpulseCenterMode[0] == 'A')
                {
                    Cfg.BCImpulseCenter = FindMaxPcm(Cfg.BCInFile2,(IFileType) Cfg.BCInFileType2[0]);
                    if (Cfg.BCImpulseCenter < 0) return 1;
                }
            ctx2.InSig = new DLReal[Cfg.BCInitWindow];
            if (ctx2.InSig == NULL) return 1;
            if (ReadSignal(Cfg.BCInFile2,ctx2.InSig,Cfg.BCInitWindow,Cfg.BCImpulseCenter,(IFileType) Cfg.BCInFileType2[0],&ctx2.PSStart,&ctx2.PSEnd) == False) return 1;
            ctx2.BCInFile = Cfg.BCInFile2;
            ctx2.PSOutFile = Cfg.PSOutFile2;
            ctx2.PSOutFileType = Cfg.PSOutFileType2;
        }
        trigger_recalculation();

	start_rest_server(8080);
	return 0;
}
