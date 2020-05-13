% Jacobian of implicit kinematic constraints of
% palh4m1IC
% projection from active to passive joints coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [8x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';
% 
% Output:
% B21 [(no of passive joints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:05
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function B21 = palh4m1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: qJ has to be [8x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_projection_jacobian_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 23:05:40
% EndTime: 2020-04-11 23:05:40
% DurationCPUTime: 0.02s
% Computational Cost: add. (103->11), mult. (130->33), div. (12->1), fcn. (98->6), ass. (0->35)
unknown=NaN(3,5);
t1 = qJ(2) + qJ(4);
t2 = cos(t1);
t3 = t2 * pkin(2);
t4 = sin(t1);
t5 = t4 * pkin(2);
t6 = sin(qJ(2));
t7 = t6 * qJ(3);
t8 = t5 - t7;
t11 = cos(qJ(2));
t12 = t11 * qJ(3);
t13 = -t3 + t12;
t16 = 0.1e1 / (-t8 * t2 * pkin(2) - t5 * t13);
t22 = cos(qJ(7));
t24 = t16 * t22 * pkin(1);
t26 = sin(qJ(7));
t28 = t16 * t26 * pkin(1);
t31 = t13 * t16;
t33 = t8 * t16;
t41 = t11 ^ 2;
t44 = t6 ^ 2;
unknown(1,1) = 0;
unknown(1,2) = (t3 * t16 * t11 + t5 * t16 * t6);
unknown(1,3) = 0;
unknown(1,4) = 0;
unknown(1,5) = (t3 * t24 + t5 * t28);
unknown(2,1) = 0;
unknown(2,2) = (t31 * t11 - t33 * t6);
unknown(2,3) = 0;
unknown(2,4) = 0;
unknown(2,5) = (t31 * t22 * pkin(1) - t33 * t26 * pkin(1));
unknown(3,1) = 0;
unknown(3,2) = (-t41 * qJ(3) * t16 - t44 * qJ(3) * t16);
unknown(3,3) = 0;
unknown(3,4) = 0;
unknown(3,5) = (-t12 * t24 - t7 * t28 + 0.1e1);
B21  = unknown;
