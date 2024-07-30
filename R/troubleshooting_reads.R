# Done in bash:

#,extract,data,from,every,fastqc,zip,file
# for,z,in,*fastqc.zip;,do,i=$(echo,"$z",|,sed,'s|fastqc.zip|fastqc/fastqc_data.txt|');echo,$i;,unzip,-p,$z,$i,>,$(echo,$i,|,tr,"/","_");done,
# 
# #,pull,all,into,a,single,file,for,plotting
# touch,N_Stats.txt
# while,read,line;,do,sed,-n,'/N-Count/{n;,:a;,/>>END_MODULE/b;,p;,n;,ba}',$line,>>,N_stats.txt;done,<,<(ls,-1,*fastqc_data.txt)

####################################



x <- read_delim("/home/gzahn/Desktop/GIT_REPOSITORIES/Dust_Data/data/raw/N_stats.txt",col_names = FALSE)

ssu_fwd <- read_delim("/home/gzahn/Desktop/GIT_REPOSITORIES/Dust_Data/data/raw/N_stats_SSU_fwd.txt",col_names = FALSE) %>% 
  mutate(amplicon="SSU",direction="FWD")
ssu_rev <- read_delim("/home/gzahn/Desktop/GIT_REPOSITORIES/Dust_Data/data/raw/N_stats_SSU_rev.txt",col_names = FALSE) %>% 
  mutate(amplicon="SSU",direction="REV")
its_fwd <- read_delim("/home/gzahn/Desktop/GIT_REPOSITORIES/Dust_Data/data/raw/N_stats_ITS_fwd.txt",col_names = FALSE) %>% 
  mutate(amplicon="ITS",direction="FWD")
its_rev <- read_delim("/home/gzahn/Desktop/GIT_REPOSITORIES/Dust_Data/data/raw/N_stats_ITS_rev.txt",col_names = FALSE) %>% 
  mutate(amplicon="ITS",direction="REV")

x <- full_join(ssu_fwd,ssu_rev) %>% 
  full_join(its_fwd) %>% 
  full_join(its_rev)
x$X1 <- factor(x$X1,
       levels = c("1","2","3","4","5","6","7","8","9","10-14","15-19","20-24",
                  "25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84",
                  "85-89","90-94","95-99","100-104","105-109","110-114","115-119","120-124","125-129","130-134","135-139","140-144",
                  "145-149","150-154","155-159","160-164","165-169","170-174","175-179","180-184","185-189","190-194","195-199","200-204",
                  "205-209","210-214","215-219","220-224","225-229","230-234","235-239","240-244","245-249","250-254","255-259","260-264",
                  "265-269","270-274","275-279","280-284","285-289","290-294","295-299","300-301"))
x %>% 
  ggplot(aes(x=X1,y=X2)) +
  geom_point(alpha=.2) +
  facet_wrap(~amplicon*direction,scales = 'free') +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5,size=8)) +
  labs(x="position",y="Count of Ns",caption = "X-position is binned.")
