fits = alpine_fit.rds subalpine_fit.rds

figs = fig/abundance-survival-time-series.pdf \
	fig/bd-load-survival.pdf \
	fig/detection-effects.pdf \
	fig/transition-plot.pdf \
	fig/mean-bd-loads.pdf \
	fig/primary-periods.pdf \
	fig/recruitment-effects.pdf \
	fig/recruitment-time-series.pdf \
	fig/survival-effects.pdf \
	fig/introduced-adult-survival.pdf
	
named_figs = figure1.pdf \
	figure2.pdf \
	figure3.pdf \
	figure4.pdf \
	figure5.pdf \
	figure6.pdf \
	figure7.pdf \
	figure8.pdf

weather_data = data/dana_meadows_snow_depth_2006-2017.csv \
	data/dana_meadows_swe_2006-2017.csv \
	data/ellery_lake_temperature_2006-2017.csv
	
appendices = appendix_s1.pdf appendix_s2.pdf

all: merged.pdf R/model-verification.pdf

merged.pdf: article.pdf $(named_figs) $(appendices)
		Rscript -e "rmarkdown::render('article.Rmd', params = list(include_figs = TRUE), output_file = 'article_w_figs.pdf')"
		pdftk article_w_figs.pdf $(appendices) cat output merged.pdf
		rm article_w_figs.pdf

article.pdf: $(named_figs) $(figs) article.Rmd library.bib
		Rscript -e "rmarkdown::render('article.Rmd')"

R/model-verification.pdf: R/model-verification.Rmd
		Rscript -e "rmarkdown::render('R/model-verification.Rmd')"

appendix_s1.pdf: appendix_s1.Rmd
		Rscript -e "rmarkdown::render('appendix_s1.Rmd')"

appendix_s2.pdf: appendix_s2.Rmd
		Rscript -e "rmarkdown::render('appendix_s2.Rmd')"



alpine_fit.rds: alpine.RData stan/uncertain-state.stan R/fit-alpine.R
		Rscript --vanilla R/fit-alpine.R

alpine.RData: R/clean-alpine-data.R data/alpine_introductions_2006-2017.csv \
	data/alpine_captures_2006-2017.csv $(weather_data)
		Rscript --vanilla R/clean-alpine-data.R

subalpine_fit.rds: subalpine.RData stan/uncertain-state.stan R/fit-subalpine.R
		Rscript --vanilla R/fit-subalpine.R

subalpine.RData: R/clean-subalpine-data.R \
	data/subalpine_introductions_2008-2017.csv \
	data/subalpine_captures_2008-2017.csv $(weather_data)
		Rscript --vanilla R/clean-subalpine-data.R
		
		

fig/abundance-survival-time-series.pdf: R/plot-abundance-time-series.R $(fits)
		Rscript --vanilla R/plot-abundance-time-series.R

fig/bd-load-survival.pdf: R/plot-bd-load-survival.R $(fits)
		Rscript --vanilla R/plot-bd-load-survival.R

fig/detection-effects.pdf: R/plot-detection-effects.R $(fits)
		Rscript --vanilla R/plot-detection-effects.R

fig/transition-plot.pdf: R/plot-transitions.R $(fits)
		Rscript --vanilla R/plot-transitions.R

fig/mean-bd-loads.pdf out/cor_prev_z.csv: R/plot-mean-bd-loads.R $(fits)
		Rscript --vanilla R/plot-mean-bd-loads.R

fig/primary-periods.pdf: R/plot-primary-periods.R $(fits)
		Rscript --vanilla R/plot-primary-periods.R

fig/recruitment-effects.pdf: R/plot-recruitment-effects.R $(fits)
		Rscript --vanilla R/plot-recruitment-effects.R

fig/recruitment-time-series.pdf: R/plot-recruitment-time-series.R $(fits)
		Rscript --vanilla R/plot-recruitment-time-series.R

fig/survival-effects.pdf: R/plot-survival-effects.R $(fits)
		Rscript --vanilla R/plot-survival-effects.R

fig/introduced-adult-survival.pdf: R/plot-introduced-adult-survival.R $(fits)
		Rscript --vanilla R/plot-introduced-adult-survival.R

figure1.pdf: fig/hmm-structure.pdf
	cp fig/hmm-structure.pdf figure1.pdf
	
figure2.pdf: fig/detection-effects.pdf
	cp fig/detection-effects.pdf figure2.pdf
	
figure3.pdf: fig/mean-bd-loads.pdf
	cp fig/mean-bd-loads.pdf figure3.pdf
	
figure4.pdf: fig/transition-plot.pdf
	cp fig/transition-plot.pdf figure4.pdf
	
figure5.pdf: fig/abundance-survival-time-series.pdf
	cp fig/abundance-survival-time-series.pdf figure5.pdf

figure6.pdf: fig/bd-load-survival.pdf
	cp fig/bd-load-survival.pdf figure6.pdf
	
figure7.pdf: fig/introduced-adult-survival.pdf
	cp fig/introduced-adult-survival.pdf figure7.pdf
	
figure8.pdf: fig/recruitment-time-series.pdf
	cp fig/recruitment-time-series.pdf figure8.pdf



clean:
		rm *.rds
		rm *.RData
		rm Rplots.pdf
		rm article.pdf
