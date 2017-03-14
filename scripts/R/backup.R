#!/home/yos2006/tools/RHEL6/R-2.15.1/bin/Rscript
###############################
library('getopt');
library('zoo');

## PARSE COMMAND LINE OPTIONS
options = getopt(
	matrix(
		c(
			'utr5', '5', 1, 'character', "5' UTR file(s)",
			'cds', 'c', 1, 'character', "CDS file(s)",
			'utr3', '3', 1, 'character', "3' UTR file(s)",
			'introns', 'i', 1, 'character', "introns file(s)",
			'flanking1kb5', 'f3', 1, 'character', "flanking 1kb upstream file(s)",
			'flanking1kb3', 'f5', 1, 'character', "flanking 1kb downstream file(s)",
			'N', 'N', 1, 'integer', "N",
			'output', 'o', 1, 'character', "Output filename",
			'legend', 'l', 1, 'character', "Legend Titles",
			'rollmean', 'r', 1, 'integer', "Rolling average distance",
			'lty', 't', 1, 'character', "R-plotting option: line type(s)",
			'col', 'col', 1, 'character', "R-plotting option: line color(s)",
			'lwd', 'w', 1, 'character', "R-plotting option: line width"
		),
		ncol = 5, 
		byrow = TRUE
	)
);

getOption = function(X, default = NULL) {
	if(!is.null(options[[X]])) {
		return(options[[X]]);
	} else {
		return(default);
	}
}

getOptions = function(X, default = NULL) {
	myopt = getOption(X);
	if(!is.null(myopt)) {
		return(unlist(strsplit(X, ",")));
	} else {
		return(default);
	}
}

N = getOption("N", 1);
output = getOption("output", "output.pdf");
colors = getOptions("col", c("blue", "red", "green4", "darkorange", "magenta", "navy"));
lwd = getOptions("lwd", 1.5);
lty = getOptions("lty", 1);
ROLLMEAN = getOption("rollmean", 5);

nrows = N + 1;
ncols = 0;
utr5_max = cds_max = introns_max = utr3_max = flanking3_max = flanking5_max = 0;

if(!is.null(options$flanking1kb5)) {
	flanking5_files = unlist(strsplit(options$flanking1kb5, ","));
	flanking5 = lapply(flanking5_files, function(X) { read.table(X); });
	flanking5_max = max(unlist(lapply(flanking5, function(X) { max(X); })));

	plotFlanking5 = TRUE;
	nBins = nrow(flanking5[[1]]);
	ncols = ncols + 1;
} else {
	plotFlanking5 = FALSE;
}

if(!is.null(options$flanking1kb3)) {
	flanking3_files = unlist(strsplit(options$flanking1kb3, ","));
	flanking3 = lapply(flanking3_files, function(X) { read.table(X); });
	flanking3_max = max(unlist(lapply(flanking3, function(X) { max(X); })));

	plotFlanking3 = TRUE;
	nBins = nrow(flanking3[[1]]);
	ncols = ncols + 1;
} else {
	plotFlanking3 = FALSE;
}

if(!is.null(options$utr5)) {
	utr5_files = unlist(strsplit(options$utr5, ","));
	utr5 = lapply(utr5_files, function(X) { read.table(X); });
	utr5_max = max(unlist(lapply(utr5, function(X) { max(X); })));
	plotUTR5 = TRUE;
	nBins = nrow(utr5[[1]]);
	ncols = ncols + 1;
} else {
	plotUTR5 = FALSE;
}

if(!is.null(options$cds)) {
	cds_files = unlist(strsplit(options$cds, ","));
	cds = lapply(cds_files, function(X) { read.table(X); });
	cds_max = max(unlist(lapply(cds, function(X) { max(X); })));
	plotCDS = TRUE;
	nBins = nrow(utr5[[1]]);
	ncols = ncols + 2 * N + 1;
} else {
	plotCDS = FALSE;
}

if(!is.null(options$introns)) {
	introns_files = unlist(strsplit(options$introns, ","));
	introns = lapply(introns_files, function(X) { read.table(X); });
	introns_max = max(unlist(lapply(introns, function(X) { max(X); })));
	plotIntrons = TRUE;
	nBins = nrow(utr5[[1]]);
	ncols = ncols + 2 * N + 1 + 2;
} else {
	plotIntrons = FALSE;
}

if(!is.null(options$utr3)) {
	utr3_files = unlist(strsplit(options$utr3, ","));
	utr3 = lapply(utr3_files, function(X) { read.table(X); });
	utr3_max = max(unlist(lapply(utr3, function(X) { max(X); })));
	plotUTR3 = TRUE;
	nBins = nrow(utr5[[1]]);
	ncols = ncols + 1 + 2;
} else {
	plotUTR3 = FALSE;
}

# TODO: FIXME
n_samples = length(utr5);
plot_max = max(utr5_max, cds_max, introns_max, utr3_max, flanking5_max, flanking3_max);
ylim = c(0, plot_max);

WIDTH = 1;
HEIGHT = 3;
W_MARGIN = 2;
H_MARGIN = 2;

plot_height = 1 / (nrows * HEIGHT + H_MARGIN * 2);
plot_width = 1 / (ncols * WIDTH + W_MARGIN * 2);
plot_margin = .1 * plot_height;

xmargin = W_MARGIN * plot_width;
ymargin = H_MARGIN * plot_height;

suffixes = c('st', 'nd', 'rd', 'th');

# PLOT

pdf(output, width = ncols * WIDTH + W_MARGIN * 2, height = nrows * HEIGHT + H_MARGIN * 2);
par(new = FALSE, mar=c(0,0,1,0), cex.main = 0.75, cex.axis = 0.75);
c_cds = c_introns = c_utr5 = c_utr3 = c_flanking5 = c_flanking3 = 1;

for(s in 1:n_samples) {
	y = 1 - plot_height * HEIGHT - ymargin - plot_margin;
	# plot boxes 0 -> N-1
	for(n in 0:(N-1)) {
		row_ncols = 0;
		if(plotFlanking5) { row_ncols = row_ncols + 1; }
		if(plotUTR5) { row_ncols = row_ncols + 1; }
		if(plotCDS) { row_ncols = row_ncols + 2 * n + 1; }
		if(plotUTR3) { row_ncols = row_ncols + 1; }
		if(plotFlanking3) { row_ncols = row_ncols + 1; }
		if(plotIntrons) { row_ncols = row_ncols + 2 * n; }
		row_xmargin = (1 - row_ncols * plot_width*WIDTH - 2*xmargin) / 2;
		row_col = colors[s %% length(colors) + 1];
		row_lty = lty[s %% length(lty) + 1];
		row_lwd = lwd[s %% length(lwd) + 1];
		
		# 5' UTR
		x = xmargin + row_xmargin;
		
		if(plotFlanking5) {
			par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
			x = x + plot_width * WIDTH;
			plot(rollmean(flanking5[[s]][,c_flanking5], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', ylim=ylim, ylab='% Peaks');
			title(paste("-1KB"));
			c_flanking5 = c_flanking5 + 1;
		}
		
		if(plotUTR5) {
			par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
			x = x + plot_width * WIDTH;
			plot(rollmean(utr5[[s]][,c_utr5], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
			title(paste("5' UTR"));
			c_utr5 = c_utr5 + 1;
		}
		
		# Plot the first N boxes
		i = 1;
		while(i <= n) {
			if(plotCDS) {
				par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
				x = x + plot_width * WIDTH;
				plot(rollmean(cds[[s]][,c_cds], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
				title(paste(i, suffixes[min(i, length(suffixes))], " CDS", sep=""));
				c_cds = c_cds + 1;
			}
			
			if(plotIntrons) {
				par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
				x = x + plot_width * WIDTH;
				plot(rollmean(introns[[s]][,c_introns], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
				title(paste(i, suffixes[min(i, length(suffixes))], " Intron", sep=""));
				c_introns = c_introns + 1;
			}
			
			i = i + 1;
		}
		
		# plot the middle CDS box
		if(plotCDS) {
			par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
			x = x + plot_width * WIDTH;
			plot(rollmean(cds[[s]][,c_cds], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
			title(paste("Middle CDS"));
			c_cds = c_cds + 1;
		}
		
		# plot the next N boxes
		i = n;
		while(i > 0) {			
			if(plotIntrons) {
				par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
				x = x + plot_width * WIDTH;
				plot(rollmean(introns[[s]][,c_introns], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
				if(i == 1) {
					title(paste("Nth Intron"));
				} else {
					title(paste("[N-", (i-1), "]th Intron", sep=""));
				}
				c_introns = c_introns + 1;
			}
			
			if(plotCDS) {
				par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
				x = x + plot_width * WIDTH;
				plot(rollmean(cds[[s]][,c_cds], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
				if(i == 1) {
					title(paste("Nth CDS"));
				} else {
					title(paste("[N-", (i-1), "]th CDS", sep=""));
				}
				c_cds = c_cds + 1;
			}
			
			i = i - 1;
		}
		
		if(plotUTR3) {
			par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
			x = x + plot_width * WIDTH;
			plot(rollmean(utr3[[s]][,c_utr3], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
			title(paste("3' UTR"));
			c_utr3 = c_utr3 + 1;
		}
		
		if(plotFlanking3) {
			par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
			x = x + plot_width * WIDTH;
			plot(rollmean(flanking3[[s]][,c_flanking3], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
			title(paste("+1KB"));
			c_flanking3 = c_flanking3 + 1;
		}
		
		y = y - plot_height * HEIGHT - 2 * plot_margin;
	}
	
	# plot row N
	n = N;
	row_ncols = 0;
	if(plotFlanking5) { row_ncols = row_ncols + 1; }
	if(plotUTR5) { row_ncols = row_ncols + 1; }
	if(plotCDS) { row_ncols = row_ncols + 2 * n + 1 + 2; }
	if(plotUTR3) { row_ncols = row_ncols + 1; }
	if(plotFlanking3) { row_ncols = row_ncols + 1; }
	if(plotIntrons) { row_ncols = row_ncols + 2 * n + 1 + 2; }
	row_xmargin = (1 - 2 * xmargin - row_ncols*plot_width*WIDTH) / 2;
	row_col = colors[s %% length(colors) + 1];
	row_lty = lty[s %% length(lty) + 1];
	row_lwd = lwd[s %% length(lwd) + 1];
	
	# 5' UTR
	x = xmargin + row_xmargin;
	
	if(plotFlanking5) {
		par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
		x = x + plot_width * WIDTH;
		plot(rollmean(flanking5[[s]][,c_flanking5], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', ylim=ylim, ylab='% Peaks');
		title(paste("-1KB"));
		c_flanking5 = c_flanking5 + 1;
	}
	
	if(plotUTR5) {
		par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
		x = x + plot_width * WIDTH;
		plot(rollmean(utr5[[s]][,c_utr5], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
		title(paste("5' UTR"));
		c_utr5 = c_utr5 + 1;
	}
	
	# Plot the first N boxes
	i = 1;
	while(i <= n) {
		if(plotCDS) {
			par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
			x = x + plot_width * WIDTH;
			plot(rollmean(cds[[s]][,c_cds], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
			title(paste(i, suffixes[min(i, length(suffixes))], " CDS", sep=""));
			c_cds = c_cds + 1;
		}
		
		if(plotIntrons) {
			par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
			x = x + plot_width * WIDTH;
			plot(rollmean(introns[[s]][,c_introns], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
			title(paste(i, suffixes[min(i, length(suffixes))], " Intron", sep=""));
			c_introns = c_introns + 1;
		}
		
		i = i + 1;
	}
	
	# plot the middle CDS box
	if(plotCDS) {
		par(fig=c(x, x+plot_width*WIDTH*3, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
		x = x + plot_width * WIDTH * 3;
		plot(rollmean(cds[[s]][,c_cds], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
		title(paste("Middle CDS"));
		c_cds = c_cds + 1;
	}
	
	# plot the middle introns box
	if(plotIntrons) {
		par(fig=c(x, x+plot_width*WIDTH*3, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
		x = x + plot_width * WIDTH * 3;
		plot(rollmean(introns[[s]][,c_introns], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
		title(paste("Middle Introns"));
		c_introns = c_introns + 1;
	}
	
	# plot the next N boxes
	i = n;
	while(i > 0) {			
		if(plotIntrons) {
			par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
			x = x + plot_width * WIDTH;
			plot(rollmean(introns[[s]][,c_introns], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
			if(i == 1) {
				title(paste("Nth Intron"));
			} else {
				title(paste("[N-", (i-1), "]th Intron", sep=""));
			}
			c_introns = c_introns + 1;
		}
		
		if(plotCDS) {
			par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
			x = x + plot_width * WIDTH;
			plot(rollmean(cds[[s]][,c_cds], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
			if(i == 1) {
				title(paste("Nth CDS"));
			} else {
				title(paste("[N-", (i-1), "]th CDS", sep=""));
			}
			c_cds = c_cds + 1;
		}
		
		i = i - 1;
	}
	
	if(plotUTR3) {
		par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
		x = x + plot_width * WIDTH;
		plot(rollmean(utr3[[s]][,c_utr3], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
		title(paste("3' UTR"));
		c_utr3 = c_utr3 + 1;
	}
	
	if(plotFlanking3) {
		par(fig=c(x, x+plot_width*WIDTH, y+plot_margin, y+plot_height*HEIGHT-plot_margin), new = TRUE);
		x = x + plot_width * WIDTH;
		plot(rollmean(flanking3[[s]][,c_flanking3], ROLLMEAN), type='l', lty=row_lty, col=row_col, lwd=row_lwd, xaxt='n', yaxt='n', ylim=ylim);
		title(paste("+1KB"));
		c_flanking3 = c_flanking3 + 1;
	}
}

#if(!is.null(options$legend)) {
#	mylegend = unlist(strsplit(options$legend, ","));
#	legend('topleft', mylegend, fill=colors[1:n_samples]);
#}
dev.off();
