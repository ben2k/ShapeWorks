#include <iostream>

// Qt includes
#include <QFileDialog>
#include <QMessageBox>

#include "SplashScreen.h"

#include <Applications/Configuration.h>
#include <Data/StudioLog.h>

#include "ui_SplashScreen.h"

namespace shapeworks {

//---------------------------------------------------------------------------
SplashScreen::SplashScreen(QWidget* parent, Preferences& preferences) :
  QDialog(parent), preferences_(preferences)
{
  this->ui_ = new Ui_SplashScreen;
  this->ui_->setupUi(this);

  this->setWindowFlags(Qt::Dialog | Qt::FramelessWindowHint);
  this->setModal(true);

  // Set up the private internals of the AppSplash class
  this->setObjectName(QString::fromUtf8("splashscreen"));

  this->ui_->version_->setText(SHAPEWORKS_VERSION);

  // Disable these since they arent being used yet.
  this->ui_->load_recent_button_->setEnabled(false);

  this->populate_recent_projects();

  connect(this->ui_->quit_button_, &QPushButton::clicked,
          this, &SplashScreen::quit);

  connect(this->ui_->recent_project_listwidget_, &QListWidget::itemPressed,
          this, &SplashScreen::enable_load_recent_button);

  connect(this->ui_->recent_project_listwidget_, &QListWidget::itemDoubleClicked,
          this, &SplashScreen::open_recent);

  connect(this->ui_->new_project_button_, &QPushButton::clicked,
          this, &SplashScreen::new_project);

  connect(this->ui_->load_recent_button_, &QPushButton::clicked,
          this, &SplashScreen::open_recent);

  connect(this->ui_->existing_project_button_, &QPushButton::clicked,
          this, &SplashScreen::open_existing);

  this->ui_->new_project_button_->setFocus();

}

//---------------------------------------------------------------------------
SplashScreen::~SplashScreen()
{
}

//---------------------------------------------------------------------------
void SplashScreen::new_project()
{
  this->close();
}

//---------------------------------------------------------------------------
void SplashScreen::quit()
{
  reinterpret_cast<QWidget*>( this->parent())->close();
}

//---------------------------------------------------------------------------
void SplashScreen::open_existing()
{

  QString filename = QFileDialog::getOpenFileName(this, tr("Open Project..."),
                                                  this->preferences_.get_last_directory(),
                                                  tr("XLSX files (*.xlsx)"));
  if (filename.isEmpty()) {
    return;
  }
  this->preferences_.set_last_directory(QFileInfo(filename).absolutePath());

  this->hide();
  QApplication::processEvents();
  emit open_project(filename);
  this->close();
}

//---------------------------------------------------------------------------
void SplashScreen::open_recent()
{

  QListWidgetItem* current_item = this->ui_->recent_project_listwidget_->currentItem();
  if (current_item == nullptr) {
    return;
  }

  QStringList user_data = current_item->data(Qt::UserRole).toStringList();

  if (user_data.size() != 2) {
    return;
  }
  QString full_file_path = user_data[0];
  QString relative_path = user_data[1];

  if (!QFileInfo::exists(full_file_path)) {
    QMessageBox::critical(0, "Project not found", "The project no longer exists.");
    return;
  }

  QDir::setCurrent(relative_path);
  this->hide();
  QApplication::processEvents();
  emit open_project(full_file_path);
  this->close();
}

//---------------------------------------------------------------------------
void SplashScreen::populate_recent_projects()
{

  QStringList recent_files = this->preferences_.get_recent_files();
  QStringList recent_paths = this->preferences_.get_recent_paths();

  int num_recent_files = qMin(recent_files.size(), (int) Preferences::MAX_RECENT_FILES);

  for (int i = 0; i < num_recent_files; i++) {
    QString text = QFileInfo(recent_files[i]).fileName();
    QListWidgetItem* new_item = new QListWidgetItem(text);
    // set the full file path and relative path as user data
    QStringList user_data = {recent_files[i], recent_paths[i]};
    new_item->setData(Qt::UserRole, QVariant(user_data));
    new_item->setToolTip(recent_files[i]);
    this->ui_->recent_project_listwidget_->addItem(new_item);
  }

}

//---------------------------------------------------------------------------
void SplashScreen::enable_load_recent_button(QListWidgetItem* not_used)
{
  this->ui_->load_recent_button_->setEnabled(true);
}

//---------------------------------------------------------------------------
void SplashScreen::resizeEvent(QResizeEvent* event)
{
  QDialog::resizeEvent(event);

  QFontMetrics fm(this->ui_->title_->font());
  int width = fm.width(this->ui_->title_->text());
  this->ui_->title_->setMinimumWidth(width);
  this->resize(width * 2, width * 1.2);
}

}

