% Jacobian of implicit kinematic constraints of
% palh4m1IC
% with respect to passive joint coordinates
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
% Phi_p [(no of constraints)x(no of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:05
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_p = palh4m1IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: qJ has to be [8x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1IC_kinconstr_impl_pas_jacobian_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobian_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 23:05:40
% EndTime: 2020-04-11 23:05:40
% DurationCPUTime: 0.01s
% Computational Cost: add. (7->4), mult. (6->4), div. (0->0), fcn. (6->4), ass. (0->16)
unknown=NaN(3,3);
t1 = qJ(2) + qJ(4);
t2 = sin(t1);
t3 = t2 * pkin(2);
t4 = sin(qJ(2));
t7 = cos(t1);
t8 = t7 * pkin(2);
t9 = cos(qJ(2));
unknown(1,1) = -t4 * qJ(3) + t3;
unknown(1,2) = t3;
unknown(1,3) = 0.0e0;
unknown(2,1) = t9 * qJ(3) - t8;
unknown(2,2) = -t8;
unknown(2,3) = 0.0e0;
unknown(3,1) = 0.1e1;
unknown(3,2) = 0.1e1;
unknown(3,3) = 0.1e1;
Phi_p  = unknown;
