"""Exercise the release comparator's failure boundaries in fresh R sessions."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "tools/replication/compare_paper_dependency.R"


@unittest.skipUnless(shutil.which("Rscript"), "Rscript unavailable")
class DependencyComparisonTests(unittest.TestCase):
    def run_r(self, code):
        return subprocess.run(
            ["Rscript", "--vanilla", "-e", code], cwd=ROOT,
            text=True, capture_output=True, check=False,
        )

    def test_quantity_budgets_and_exact_discrete_decisions(self):
        result = self.run_r(f'''
          source({str(SCRIPT)!r})
          budget <- list(atol=1e-9, rtol=1e-12, scale=1)
          stopifnot(compare_numeric(c(0,1), c(1e-10,1), budget)$pass)
          stopifnot(!compare_numeric(0, 1e-5, budget)$pass)
          stopifnot(!compare_numeric(c(1,NA), c(1,0), budget)$pass)
          stopifnot(!compare_discrete(c(1L,2L), c(2L,1L))$pass)
          stopifnot(compare_discrete(c(1L,2L), c(1L,2L))$pass)
          stopifnot(!compare_numeric(0, -0, budget)$bitwise)
        ''')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_identity_cannot_substitute_a_traversal_path_for_the_loaded_library(self):
        result = self.run_r(f'''
          source({str(SCRIPT)!r})
          pkg <- tryCatch(find.package("RdistNN"),error=function(e)NULL)
          if (is.null(pkg)) quit(status=77)
          file <- tempfile(fileext=".json")
          jsonlite::write_json(list(version=read.dcf(file.path(pkg,"DESCRIPTION"))[1,"Version"],
            source_identity="invalid native identity control",
            files=setNames(list(file_sha256(file.path(pkg,"DESCRIPTION"))),"libs/../DESCRIPTION")),
            file,auto_unbox=TRUE)
          err <- tryCatch(load_identified_package(dirname(pkg),file),error=function(e)e)
          stopifnot(inherits(err,"error"),grepl("canonical relative",conditionMessage(err)))
        ''')
        if result.returncode == 77:
            self.skipTest("Installed RdistNN unavailable")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_refuses_existing_output_and_repeated_arguments(self):
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run(
                ["Rscript", "--vanilla", str(SCRIPT),
                 f"--output-dir={directory}", f"--output-dir={directory}"],
                text=True, capture_output=True, check=False,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Repeated argument", result.stderr)
            result = self.run_r(f'''
              source({str(SCRIPT)!r})
              new_output({directory!r})
            ''')
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("already exists", result.stderr)

    def test_frozen_packet_tampering_is_rejected(self):
        result = self.run_r(f'''
          source({str(SCRIPT)!r})
          path <- tempfile(); writeLines("changed", path)
          expect <- paste(rep("0",64),collapse="")
          err <- tryCatch(checked_hash(path, expect), error=identity)
          stopifnot(inherits(err,"error"), grepl("Hash mismatch",conditionMessage(err)))
        ''')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_same_installation_control_uses_distinct_processes_and_verified_engines(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_r(f'''
              source({str(SCRIPT)!r})
              if (!requireNamespace("RdistNN",quietly=TRUE)) quit(status=77)
              root <- {directory!r}
              pkg <- find.package("RdistNN")
              files <- list.files(pkg,recursive=TRUE,full.names=TRUE)
              files <- files[!dir.exists(files)]
              hashes <- as.list(setNames(vapply(files,file_sha256,character(1)),
                substring(files,nchar(pkg)+2L)))
              identity <- list(version=as.character(packageVersion("RdistNN")),
                source_identity="same-installation validation control",files=hashes)
              jsonlite::write_json(identity,file.path(root,"identity.json"),auto_unbox=TRUE)
              x <- list(id="paired",X=matrix(1:4,ncol=1),X0=matrix(0),
                PHI=cbind(c(2,1,4,3),c(1,-1,2,1),-c(2,1,4,3)),
                ranks=matrix(1:4),wn=c(1/2,1/3,1/6,0),wnm1=c(2/3,1/3,0),
                row_ids=1:4,target_ids="zero",k=1L,s=2L,geometry="l2",display_digits=6L)
              saveRDS(x,file.path(root,"packet.rds"))
              fields <- c("full_weights","deletion_weights","full_mass","deletion_mass","ess",
                "components","deleted","deleted_averages","centered","pseudo","component_covariance",
                "cate_target_covariance","cate_estimate","cate_se","estimate","se","lower","upper",
                "local_limit_variance","distances","wrapper_cate_estimate","wrapper_cate_se",
                "wrapper_lower","wrapper_upper")
              budgets <- setNames(lapply(fields,function(f) list(atol=1e-12,rtol=1e-12,scale=1)),fields)
              manifest <- list(schema_version=1L,paper_files=list(),cases=list(list(id=x$id,
                file="packet.rds",sha256=file_sha256(file.path(root,"packet.rds")),
                scope="validation control",budgets=budgets)))
              jsonlite::write_json(manifest,file.path(root,"cases.json"),auto_unbox=TRUE)
              args <- c("--vanilla",{str(SCRIPT)!r},paste0("--paper=",root),
                paste0("--historical-library=",dirname(pkg)),paste0("--candidate-library=",dirname(pkg)),
                paste0("--historical-identity=",root,"/identity.json"),
                paste0("--candidate-identity=",root,"/identity.json"),
                paste0("--cases=",root,"/cases.json"),paste0("--output-dir=",root,"/result"))
              status <- system2(file.path(R.home("bin"),"Rscript"),shQuote(args))
              if (status != 0) {{
                cat(readLines(file.path(root,"result/historical/process.log")),sep="\\n")
                stop("Self-control failed")
              }}
              a <- readRDS(file.path(root,"result/historical/result.rds"))
              b <- readRDS(file.path(root,"result/candidate/result.rds"))
              stopifnot(a$identity$pid != b$identity$pid,
                a$identity$engine_calls$compiled > 0L,a$identity$engine_calls$forced_R == 0L,
                identical(a$modes,b$modes))
              stopifnot(identical(names(a$modes$compiled$paired),
                c("A_arithmetic","B_exact_weights","C_geometry","D_wrapper")),
                identical(a$modes$compiled$paired$D_wrapper$ranks,matrix(1:4)))
              err <- tryCatch(load_identified_package(dirname(pkg),file.path(root,"identity.json")),error=function(e)e)
              stopifnot(inherits(err,"error"),grepl("already loaded",conditionMessage(err)))
            ''')
            if result.returncode == 77:
                self.skipTest("Installed RdistNN unavailable for fresh-process control")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
