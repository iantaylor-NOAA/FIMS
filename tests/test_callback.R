
library(FIMS)


recruitment_call<-function(rzero,phi_0,steep, spawners,ssbzero){

    recruits<-(0.8 * rzero * steep * spawners) /
                  (0.2 * phi_0 * rzero * (1.0 - steep) + spawners * (steep - 0.2));
    return(recruits);
}


# Recruitment
recruitment <- new(RecruitmentRCallback)
p1 <- new(Parameter)
p1$value = 3;
p2 <- new(Parameter)
p2$value = 4;
p3 <- new(Parameter)
p3$value = 4;
recruitment$AddArgument("rzero", p1)
recruitment$AddArgument("phi_0", p2)
recruitment$AddArgument("steep", p2)
recruitment$SetFunction(recruitment_call)
recruitment$evaluate(10.0,2.0)
