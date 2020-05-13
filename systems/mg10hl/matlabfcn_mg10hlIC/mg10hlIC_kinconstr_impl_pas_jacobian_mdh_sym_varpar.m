% Jacobian of implicit kinematic constraints of
% mg10hlIC
% with respect to passive joint coordinates
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
% Phi_p [(no of constraints)x(no of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 13:09
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_p = mg10hlIC_kinconstr_impl_pas_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'mg10hlIC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'mg10hlIC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobian_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 13:09:18
% EndTime: 2020-04-11 13:09:18
% DurationCPUTime: 0.03s
% Computational Cost: add. (13->8), mult. (12->10), div. (0->0), fcn. (12->10), ass. (0->63)
unknown=NaN(7,7);
t1 = qJ(2) + pkin(14) + qJ(3);
t2 = cos(t1);
t4 = sin(qJ(10));
t6 = sin(t1);
t8 = cos(qJ(10));
t10 = sin(qJ(5));
t12 = qJ(9) + qJ(11);
t13 = sin(t12);
t14 = t13 * qJ(13);
t15 = sin(qJ(9));
t18 = cos(qJ(5));
t20 = cos(t12);
t21 = t20 * qJ(13);
t22 = cos(qJ(9));
unknown(1,1) = t2 * pkin(3);
unknown(1,2) = 0.0e0;
unknown(1,3) = 0.0e0;
unknown(1,4) = 0.0e0;
unknown(1,5) = -t4 * pkin(4);
unknown(1,6) = 0.0e0;
unknown(1,7) = 0.0e0;
unknown(2,1) = t6 * pkin(3);
unknown(2,2) = 0.0e0;
unknown(2,3) = 0.0e0;
unknown(2,4) = 0.0e0;
unknown(2,5) = t8 * pkin(4);
unknown(2,6) = 0.0e0;
unknown(2,7) = 0.0e0;
unknown(3,1) = 0.0e0;
unknown(3,2) = 0.0e0;
unknown(3,3) = -t10 * pkin(7);
unknown(3,4) = t15 * pkin(5) - t14;
unknown(3,5) = 0.0e0;
unknown(3,6) = -t14;
unknown(3,7) = 0.0e0;
unknown(4,1) = 0.0e0;
unknown(4,2) = 0.0e0;
unknown(4,3) = t18 * pkin(7);
unknown(4,4) = -t22 * pkin(5) + t21;
unknown(4,5) = 0.0e0;
unknown(4,6) = t21;
unknown(4,7) = 0.0e0;
unknown(5,1) = -0.1e1;
unknown(5,2) = 0.0e0;
unknown(5,3) = 0.0e0;
unknown(5,4) = 0.0e0;
unknown(5,5) = 0.1e1;
unknown(5,6) = 0.0e0;
unknown(5,7) = 0.1e1;
unknown(6,1) = 0.0e0;
unknown(6,2) = 0.0e0;
unknown(6,3) = -0.1e1;
unknown(6,4) = 0.1e1;
unknown(6,5) = 0.0e0;
unknown(6,6) = 0.1e1;
unknown(6,7) = 0.0e0;
unknown(7,1) = 0.0e0;
unknown(7,2) = -0.1e1;
unknown(7,3) = 0.0e0;
unknown(7,4) = -0.1e1;
unknown(7,5) = 0.0e0;
unknown(7,6) = 0.0e0;
unknown(7,7) = 0.0e0;
Phi_p  = unknown;
