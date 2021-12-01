#pragma once

#include <string>

#include <QWidget>
#include <QSharedPointer>
#include <QProgressDialog>
#include <QElapsedTimer>
#include <QObject>

#include <Data/Preferences.h>

class Ui_DataTool;

namespace shapeworks {

class Session;
class ShapeWorksStudioApp;

class DataTool : public QWidget {
  Q_OBJECT;
public:

  DataTool(Preferences& prefs);
  ~DataTool();

  //! Set the pointer to the session
  void set_session(QSharedPointer<Session> session);

  //! activate this tool
  void activate();

  void disable_actions();

  void enable_actions();

  void update_table();

  void update_notes();

  std::string get_notes();

public Q_SLOTS:

  //void on_table_open_button_toggled();

  void delete_button_clicked();

Q_SIGNALS:
  void import_button_clicked();

private:

  Preferences& preferences_;

  Ui_DataTool* ui_;
  QSharedPointer<Session> session_;
};
}
