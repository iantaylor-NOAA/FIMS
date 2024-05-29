
test_that("Parameter vector works as expected", {

v_size <- 10
v1_value <- 1.0
v2_value <- 2.0

#Test that default constructor works
v0 <- new(ParameterVector)
expect_equal(length(v0), 1)
expect_equal(v0$at(1)$value, 0)

# Test that constructor that initializes based on size works.
v1 <- new(ParameterVector, v_size)
v1$fill(v1_value)
for(i in 1:v_size){
   expect_equal(v1$at(i)$value, v1_value)
}

# Test that constructor that takes vector and size works.
v2 <- new(ParameterVector, rep(v2_value, v_size), v_size)
for(i in 1:v_size){
   expect_equal(v2$at(i)$value, v2_value)
}

#Test that plus operator works
v3 <- v1+v2
v3_value<-v1_value+v2_value
for(i in 1:v_size){
   expect_equal(v3$at(i)$value, v3_value)
}

#Test that minus operator works
v3 <- v1-v2
v3_value <- v1_value-v2_value

for(i in 1:v_size){
   expect_equal(v3$at(i)$value, v3_value)
}

#Test that multiply operator works
v3 <- v1*v2
v3_value <- v1_value*v2_value

for(i in 1:v_size){
   expect_equal(v3$at(i)$value, v3_value)
}

#Test that divide operator works
v3<- v1/v2
v3_value<-v1_value/v2_value

for(i in 1:v_size){
   expect_equal(v3$at(i)$value, v3_value)
}

v_acos_test <- acos(v1)
for(i in 1:v_size){
   expect_equal(v_acos_test$at(i)$value, acos(v1_value))
}


v_asin_test <- asin(v1)
for(i in 1:v_size){
   expect_equal(v_asin_test$at(i)$value, asin(v1_value))
}


v_atan_test <- atan(v1)
for(i in 1:v_size){
   expect_equal(v_atan_test$at(i)$value, atan(v1_value))
}


v_cos_test <- cos(v1)
for(i in 1:v_size){
   expect_equal(v_cos_test$at(i)$value, cos(v1_value))
}


v_cosh_test <- cosh(v1)
for(i in 1:v_size){
   expect_equal(v_cosh_test$at(i)$value, cosh(v1_value))
}


v_sin_test <- sin(v1)
for(i in 1:v_size){
   expect_equal(v_sin_test$at(i)$value, sin(v1_value))
}


v_sinh_test <- sinh(v1)
for(i in 1:v_size){
   expect_equal(v_sinh_test$at(i)$value, sinh(v1_value))
}


v_tan_test <- tan(v1)
for(i in 1:v_size){
   expect_equal(v_tan_test$at(i)$value, tan(v1_value))
}


v_tanh_test <- tanh(v1)
for(i in 1:v_size){
   expect_equal(v_tanh_test$at(i)$value, tanh(v1_value))
}


v_exp_test <- exp(v1)
for(i in 1:v_size){
   expect_equal(v_exp_test$at(i)$value, exp(v1_value))
}


v_log10_test <- log10(v1)
for(i in 1:v_size){
   expect_equal(v_log10_test$at(i)$value, log10(v1_value))
}


v_sqrt_test <- sqrt(v1)
for(i in 1:v_size){
   expect_equal(v_sqrt_test$at(i)$value, sqrt(v1_value))
}


v_log_test <- log(v1)
for(i in 1:v_size){
   expect_equal(v_log_test$at(i)$value, log(v1_value))
}





#Test that created IDs are unique
p<-new(ParameterVector, 100)

p[1]$value
p[1]$value <- 1
p[1]$value
k <- p[1]$id
for(i in 2:length(p)){
   expect_equal(p[i]$id, k + 1)
   k <- p[i]$id
}

#Test that resize works
p$resize(5)
for( i in 1:length(p)){
   print(p[i]$id)
}

p$resize(10)

for( i in 1:length(p)){
print(p[i]$id)
}


var<-p[1] + cos(p[2])
#str(var)
print(p[1]$value)
print(p[2]$value)
print(cos(p[2]$value))
#str(var)
#par$value<-3.1459
#str(var)
#a<-apply(X = v, MARGIN = 1, FUN = sum)
cc<-c(1,2,3)
dim(cc)
l<-p$data
r<-sum(p)
#str(r)
#q()
#str(l)
#a<-lapply(X = l, MARGIN = 1, FUN = sum)
x <- vector("numeric",   length = 10)


typeof(x)
})
