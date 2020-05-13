% Jacobian of implicit kinematic constraints of
% mg10hlIC
% with respect to active joint coordinates
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
% Phi_a [(no of constraints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 13:09
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_a = mg10hlIC_kinconstr_impl_act_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'mg10hlIC_kinconstr_impl_act_jacobian_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'mg10hlIC_kinconstr_impl_act_jacobian_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobian_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 13:09:18
% EndTime: 2020-04-11 13:09:18
% DurationCPUTime: 0.02s
% Computational Cost: add. (10->6), mult. (4->4), div. (0->0), fcn. (6->6), ass. (0->51)
unknown=NaN(7,6);
t1 = qJ(2) + pkin(14) + qJ(3);
t2 = cos(t1);
t4 = qJ(2) + pkin(14);
t5 = cos(t4);
t8 = sin(t1);
t10 = sin(t4);
t13 = qJ(9) + qJ(11);
t14 = cos(t13);
t15 = sin(t13);
unknown(1,1) = 0;
unknown(1,2) = (pkin(1) * t5 + t2 * pkin(3));
unknown(1,3) = 0;
unknown(1,4) = 0;
unknown(1,5) = 0;
unknown(1,6) = 0;
unknown(2,1) = 0;
unknown(2,2) = (pkin(1) * t10 + t8 * pkin(3));
unknown(2,3) = 0;
unknown(2,4) = 0;
unknown(2,5) = 0;
unknown(2,6) = 0;
unknown(3,1) = 0;
unknown(3,2) = 0;
unknown(3,3) = 0;
unknown(3,4) = 0;
unknown(3,5) = 0;
unknown(3,6) = t14;
unknown(4,1) = 0;
unknown(4,2) = 0;
unknown(4,3) = 0;
unknown(4,4) = 0;
unknown(4,5) = 0;
unknown(4,6) = t15;
unknown(5,1) = 0;
unknown(5,2) = -1;
unknown(5,3) = 0;
unknown(5,4) = 0;
unknown(5,5) = 0;
unknown(5,6) = 0;
unknown(6,1) = 0;
unknown(6,2) = 0;
unknown(6,3) = 0;
unknown(6,4) = 0;
unknown(6,5) = 0;
unknown(6,6) = 0;
unknown(7,1) = 0;
unknown(7,2) = 0;
unknown(7,3) = 0;
unknown(7,4) = 0;
unknown(7,5) = 0;
unknown(7,6) = 0;
Phi_a  = unknown;
