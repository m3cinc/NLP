# spell3.R
spell_more<-list()
spell<-c("one","two","three","four","five","six","seven","eight","nine","ten",
         "eleven","twelve","thirteen","fourteen","fifteen","sixteen","seventeen",
         "eighteen","nineteen")
tens<-c("twenty","thirty","forty","fifty","sixty","seventy","eighty","ninety") 
for (j in 1:8) {
        spell<-c(spell,tens[j])
        spell_more<-sapply(1:9,function(i){ spell_more[i]<-spell[i]})
        spell_more<-paste0(tens[j],"-",spell_more)
        spell<-c(spell,spell_more)
}
hundreds<-list()
for (j in 1:9) {
        hundreds<-c(hundreds,paste0(spell[j]," hundred "))
        spell<-c(spell,hundreds[j])
        for (l in 1:19) {
                spell<-c(spell,paste0(hundreds[j],spell[l]))
        }
        for (l in 1:8) {
        spell<-c(spell,paste0(hundreds[j],tens[l]))
        spell_more<-sapply(1:9,function(i){ spell_more[i]<-spell[i]})
        spell_more<-paste0(hundreds[j],tens[l],"-",spell_more)
        spell<-c(spell,spell_more)
        }
}
spell <-c ("zero", spell, "one thousand" )
spell_number <- array(0L:1000L)
spell_number <- setNames(spell_number,spell)
saveRDS(spell_number,file = "spell_number.RDS")
rm(spell_more,spell,tens,hundreds)
