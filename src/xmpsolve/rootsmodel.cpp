#include "rootsmodel.h"
#include <QDebug>

namespace xmpsolve {

RootsModel::RootsModel(QObject *parent) :
    QAbstractListModel(parent)
{
    length = 0;
}

QHash<int, QByteArray>
RootsModel::roleNames() const
{
    QHash<int, QByteArray> role_names;

    role_names.insert(RADIUS, "radius");
    role_names.insert(STATUS, "status");
    role_names.insert(SHORT_APPROXIMATION, "short_approximation");
    role_names.insert(ROOT, "root");

    return role_names;
}

int
RootsModel::rowCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent);
    return length;
}

QVariant
RootsModel::data(const QModelIndex &index, int role) const
{
    int i = index.row();

    if (role != Qt::DisplayRole && role < Qt::UserRole)
        return QVariant();

    if (i < 0 || i > length)
        return QVariant();
    else
    {
        rdpe_t root_module;
        mpc_rmod (root_module, m_roots.at(i)->value);
        int digits = (rdpe_Esp (root_module) - rdpe_Esp (m_roots.at(i)->radius)) / LOG2_10;
        char * buffer = NULL;

        switch (role)
        {
            case SHORT_APPROXIMATION:
                digits = 4;
                // fallthrough
            case Qt::DisplayRole:
                buffer = new char[2 * digits + 15];
                if (m_roots[i]->get_imag_part() > 0)
                    gmp_sprintf (buffer, "%.*Ff + %.*Ffi", digits, mpc_Re (m_roots[i]->value),
                                 digits, mpc_Im (m_roots[i]->value));
                else
                    gmp_sprintf (buffer, "%.*Ff %.*Ffi", digits, mpc_Re (m_roots[i]->value),
                                 digits, mpc_Im (m_roots[i]->value));
                return QString(buffer);
                break;

            case STATUS:
                return QString(MPS_ROOT_STATUS_TO_STRING (m_roots[i]->status));

            case RADIUS:
                return QString("%1").arg(m_roots[i]->get_radius());

            case ROOT:
                return QVariant::fromValue((void*) m_roots[i]);

            default:
                qDebug() << "Invalid role";
                return QVariant();
        }
    }
}

void
RootsModel::setRoots(QList<Root *> roots)
{
    beginResetModel();

    length = 0;
    m_roots = roots;
    length = roots.length();

    endResetModel();
}



} // namespace xmpsolve
