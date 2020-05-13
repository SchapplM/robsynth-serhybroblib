% Implicit kinematic constraints of
% mg10hlIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% h [(no of constraints)x1]
%   Implicit constraint equations (e.g. closed loops)
%   In a valid robot configuration qJ this has to be zero

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 13:09
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function h = mg10hlIC_kinconstr_impl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'mg10hlIC_kinconstr_impl_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'mg10hlIC_kinconstr_impl_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 13:09:18
% EndTime: 2020-04-11 13:09:18
% DurationCPUTime: 0.04s
% Computational Cost: add. (27->23), mult. (13->13), div. (0->0), fcn. (12->12), ass. (0->22)
unknown=NaN(7,1);
t1 = cos(qJ(10));
t3 = qJ(2) + pkin(14) + qJ(3);
t4 = sin(t3);
t6 = qJ(2) + pkin(14);
t7 = sin(t6);
t10 = sin(qJ(10));
t12 = cos(t3);
t14 = cos(t6);
t17 = cos(qJ(5));
t19 = qJ(9) + qJ(11);
t20 = cos(t19);
t22 = cos(qJ(9));
t25 = sin(qJ(5));
t27 = sin(t19);
t29 = sin(qJ(9));
unknown(1,1) = pkin(1) * t7 + t4 * pkin(3) + t1 * pkin(4) - pkin(12) - pkin(13);
unknown(2,1) = -pkin(1) * t14 - t12 * pkin(3) + t10 * pkin(4);
unknown(3,1) = -t22 * pkin(5) + t17 * pkin(7) + t20 * qJ(13) + pkin(6);
unknown(4,1) = -t29 * pkin(5) + t25 * pkin(7) + t27 * qJ(13);
unknown(5,1) = -pi / 0.2e1 + qJ(10) + qJ(12) - qJ(2) - pkin(14) - qJ(3);
unknown(6,1) = -qJ(5) + qJ(9) + qJ(11);
unknown(7,1) = -qJ(4) - qJ(9);
h  = unknown;
