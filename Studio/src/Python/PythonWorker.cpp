// pybind
#include <pybind11/embed.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using namespace pybind11::literals; // to bring in the `_a` literal


#include <sstream>
#include <iostream>

#include <QProcess>
#include <QMessageBox>

#include <Groom/QGroom.h>
#include <Data/Shape.h>
#include <Python/PythonWorker.h>
#include <DeepSSM/QDeepSSM.h>
#include <Data/StudioLog.h>
#include <Libs/Optimize/Optimize.h>
#include <Libs/Optimize/OptimizeParameters.h>

namespace shapeworks {

//---------------------------------------------------------------------------
class PythonLogger {
public:

  void set_callback(const std::function<void(std::string)>& callback)
  {
    this->callback_ = callback;
  }

  void set_progress_callback(const std::function<void(double)>& callback)
  {
    this->progress_callback_ = callback;
  }

  void cpp_log(std::string msg)
  {
    this->callback_(msg);
  }

  void cpp_progress(double progress)
  {
    this->progress_callback_(progress * 100);
  }

  void clear_abort()
  {
    this->aborted_ = false;
  }

  void set_abort()
  {
    this->aborted_ = true;
  }

  bool check_abort()
  {
    return this->aborted_;
  }

private:
  std::function<void(std::string)> callback_;
  std::function<void(double)> progress_callback_;
  std::atomic<bool> aborted_{false};
};

//---------------------------------------------------------------------------
PYBIND11_EMBEDDED_MODULE(logger, m)
{
  py::class_<PythonLogger, std::shared_ptr<PythonLogger>>(m, "PythonLogger")
    .def(py::init<>())
    .def("log", &PythonLogger::cpp_log)
    .def("check_abort", &PythonLogger::check_abort)
    .def("progress", &PythonLogger::cpp_progress);
};

//---------------------------------------------------------------------------
PythonWorker::PythonWorker()
{
  this->thread_ = new QThread(this);
  this->moveToThread(this->thread_);
  connect(thread_, &QThread::started, this, &PythonWorker::init, Qt::QueuedConnection);
  this->thread_->start();
}

//---------------------------------------------------------------------------
PythonWorker::~PythonWorker()
{
  this->end_python();
  this->thread_->wait();
  delete this->thread_;
}

//---------------------------------------------------------------------------
void PythonWorker::set_deep_ssm(QSharedPointer<QDeepSSM> deep_ssm)
{
  this->deep_ssm_ = deep_ssm;
  this->deep_ssm_->moveToThread(this->thread_);
}

//---------------------------------------------------------------------------
void PythonWorker::start_deepssm_augmentation()
{
  this->deep_ssm_->run_augmentation();
  this->finish_job();
}

//---------------------------------------------------------------------------
void PythonWorker::start_deepssm_training()
{
  this->deep_ssm_->run_training();
  this->finish_job();
}

//---------------------------------------------------------------------------
void PythonWorker::run_job(PythonWorker::JobType job)
{
  if (job == PythonWorker::JobType::DeepSSM_AugmentationType) {
    QMetaObject::invokeMethod(this, "start_deepssm_augmentation");
  }
  else if (job == PythonWorker::JobType::DeepSSM_TrainingType) {
    QMetaObject::invokeMethod(this, "start_deepssm_training");
  }
}

//---------------------------------------------------------------------------
void PythonWorker::init()
{
  std::string home = getenv("HOME");
#ifdef _WIN32
  home = getenv("USERPROFILE");
#endif

  // read list generated by something like this:
  // python -c "import sys; print('\n'.join(sys.path))" > $HOME/.shapeworks/python_path.txt
  std::vector<std::string> python_path;
  std::fstream file;
  file.open(home + "/.shapeworks/python_path.txt", std::ios::in);
  if (file.is_open()) {
    std::string line;
    while (getline(file, line)) {
      python_path.push_back(line);
    }
    file.close();
  }
  else {
#ifdef _WIN32
    emit error_message(QString::fromStdString("Unable to initialize Python\nPlease run install_shapeworks.bat");
#else
    emit error_message(
      QString::fromStdString("Unable to initialize Python\nPlease run install_shapeworks.sh"));
#endif
  }

  try {
    py::initialize_interpreter();

    py::module sys = py::module::import("sys");

#ifdef __APPLE__
    setenv("OMP_NUM_THREADS", "1", 1);
#endif

    sys.attr("path") = python_path;

    // this is necessary or the plots will crash the process
    py::module py_matplot_lib = py::module::import("matplotlib");
    py_matplot_lib.attr("use")("agg");

    this->python_logger_ = QSharedPointer<PythonLogger>::create();
    this->python_logger_->set_callback(
      std::bind(&PythonWorker::incoming_python_message, this, std::placeholders::_1));
    this->python_logger_->set_progress_callback(
      std::bind(&PythonWorker::incoming_python_progress, this, std::placeholders::_1));
    py::module logger = py::module::import("logger");

    py::module sw_utils = py::module::import("shapeworks.utils");
    py::object set_sw_logger = sw_utils.attr("set_sw_logger");
    set_sw_logger(this->python_logger_.data());

    STUDIO_LOG_MESSAGE("Embedded Python Interpreter Initialized");
  } catch (py::error_already_set& e) {
    emit error_message(QString::fromStdString("Error initializing Python:\n") + e.what());
  } catch (const std::exception& e) {
    emit error_message(QString::fromStdString("Error initializing Python:\n") + e.what());
  }
}

//---------------------------------------------------------------------------
void PythonWorker::incoming_python_message(std::string message_string)
{
  emit message(QString::fromStdString(message_string));
}

//---------------------------------------------------------------------------
void PythonWorker::end_python()
{
  QMetaObject::invokeMethod(this, "finalize_python");
}

//---------------------------------------------------------------------------
void PythonWorker::finalize_python()
{
  py::finalize_interpreter();
  this->thread_->exit();
}

//---------------------------------------------------------------------------
void PythonWorker::incoming_python_progress(double value)
{
  emit progress(value);
}

//---------------------------------------------------------------------------
void PythonWorker::abort_job()
{
  this->python_logger_->set_abort();
}

//---------------------------------------------------------------------------
void PythonWorker::finish_job()
{
  this->python_logger_->clear_abort();
  emit job_finished();
}

}
