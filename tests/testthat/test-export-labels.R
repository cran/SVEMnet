label_worklist <- function() {
  data.frame(recipe_id=c("r2","r1"),slot=c(2L,1L),role=c("optimum","incumbent_repeat"),
    nominated_by=c("forward","observed"),x=c(.2,.1),y_pred=c(5,4))
}

test_that("explicit export labels distinguish preparations from physical recipes", {
  work <- label_worklist()
  original <- work
  plain <- svem_export_candidates_csv(work,write_file=FALSE)
  labelled <- svem_export_candidates_csv(work,write_file=FALSE,
    planned_date="2026-11-02",preparation_prefix="studyA-20261102-round1")
  expect_identical(work,original)
  expect_identical(labelled$recipe_id,plain$recipe_id)
  expect_identical(labelled[names(plain)],plain)
  expect_identical(labelled$preparation_id,c("studyA-20261102-round1-001","studyA-20261102-round1-002"))
  expect_equal(labelled$planned_date,rep(as.Date("2026-11-02"),2))
  expect_identical(labelled,svem_export_candidates_csv(work,write_file=FALSE,
    planned_date=as.Date("2026-11-02"),preparation_prefix="studyA-20261102-round1"))
  later <- svem_export_candidates_csv(work,write_file=FALSE,
    planned_date="2026-11-02",preparation_prefix="studyA-20261102-round2")
  expect_false(any(later$preparation_id %in% labelled$preparation_id))
  expect_identical(later$recipe_id,labelled$recipe_id)
  expect_error(svem_export_candidates_csv(labelled,write_file=FALSE,
    planned_date="2026-11-03"),"already has planned_date")
  expect_error(svem_export_candidates_csv(labelled,write_file=FALSE,
    preparation_prefix="new"),"already has preparation_id")
  # User-defined preparation IDs survive worklist ordering.
  work$preparation_id <- c("custom-b","custom-a")
  expect_identical(svem_export_candidates_csv(work,write_file=FALSE)$preparation_id,
                   c("custom-a","custom-b"))
})

test_that("optional labels preserve empty exports and validate explicit date/prefix", {
  empty <- label_worklist()[FALSE,]
  out <- svem_export_candidates_csv(empty,write_file=FALSE,
    planned_date="2026-11-02",preparation_prefix="test-batch")
  expect_equal(nrow(out),0)
  expect_true(all(c("planned_date","preparation_id") %in% names(out)))
  expect_s3_class(out$planned_date,"Date")
  expect_identical(out$preparation_id,character())
  for (bad in list("02/11/2026","26-11-02","2026-02-30","2026-11-02 extra",NA_character_,c("2026-11-01","2026-11-02"),1)) {
    expect_error(svem_export_candidates_csv(label_worklist(),write_file=FALSE,
      planned_date=bad),"planned_date")
  }
  for (bad in list("",NA_character_,"a\nb",1,c("a","b"))) {
    expect_error(svem_export_candidates_csv(label_worklist(),write_file=FALSE,
      preparation_prefix=bad),"preparation_prefix")
  }
})

test_that("legacy selections change only when labeling is explicitly requested", {
  tab <- data.frame(x=1:8,score=seq(.1,.8,length.out=8))
  selection <- svem_select_from_score_table(tab,target="score",k=0)
  plain <- svem_export_candidates_csv(selection,write_file=FALSE)
  labelled <- svem_export_candidates_csv(selection,write_file=FALSE,
    planned_date="2026-11-02",preparation_prefix="legacy-round1")
  expect_identical(labelled[names(plain)],plain)
  expect_identical(names(labelled),c(names(plain),"planned_date","preparation_id"))
  expect_identical(labelled$preparation_id,"legacy-round1-001")
  expect_equal(labelled$planned_date,as.Date("2026-11-02"))
  expect_identical(svem_export_candidates_csv(selection,write_file=FALSE,
    planned_date=NULL,preparation_prefix=NULL),plain)
})

